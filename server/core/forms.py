from django import forms
from django.contrib.auth.password_validation import password_validators_help_text_html
from django.core.validators import ValidationError
from django.contrib.auth.password_validation import validate_password
from django.contrib.auth.models import User
import warnings

import re
import io
from Bio import SeqIO, PDB
from Bio.PDB.PDBExceptions import PDBConstructionWarning


def validate_sequence(input):
    aa1 = 'A R N D C Q E G H I L K M F P S T W Y V X'
    aa1_split = aa1.split()

    for x in input:
        if x not in aa1_split:
            raise ValidationError('Protein sequence must belong to 21 letter alphabet (found "%s"):\n%s' % (x, aa1))


def validate_fasta(input):
    fasta = io.StringIO(input)
    try:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
    except Exception as e:
        raise ValidationError('Wrong format of FASTA input')
    fasta.close()

    if len(fasta_dict) > 1:
        raise ValidationError('More, than one sequence provided')

    if len(fasta_dict) < 1:
        raise ValidationError('No sequences provided')

    return fasta_dict


def validate_file(f):
    allowed_size_mb = 8
    if f.size > 1024 * 1024 * allowed_size_mb:
        raise ValidationError(f'Uploaded file cannot be larger than {allowed_size_mb} MB')


def validate_pdb(pdb_file):
    if not pdb_file.name.endswith('.pdb'):
        raise ValidationError('Must be a file with .pdb extension')

    warnings.simplefilter('ignore', PDBConstructionWarning)  # Silence annoying warnings cluttering web server's log
    parser = PDB.PDBParser()

    # if no decoding is done, parser just skips it all
    with io.StringIO(pdb_file.read().decode("utf-8")) as buf:
        pdb_file.seek(0)

        try:
            pdb = parser.get_structure('pdb', buf)
        except Exception as e:
            raise ValidationError("Couldn't parse PDB file: " + str(e))

        if len(list(pdb.get_models())) > 1:
            raise ValidationError("Only PDB files containing one model are supported")

        if len(list(pdb.get_atoms())) == 0:
            raise ValidationError("Didn't find any atoms in PDB file")

        for chain in pdb.get_chains():
            chain_id = chain.get_id()[-1]
            if len(chain_id) != 1:
                raise ValidationError("Invalid chain is present: " + str(chain_id))
            if chain_id == ' ':
                raise ValidationError("PDB can not contain atoms with empty chain field")

        for residue in pdb.get_residues():
            hetero_flag, sequence_num, insertion_code = residue.get_id()
            if insertion_code != ' ':
                res_id = f"{residue.get_resname()} {sequence_num}{insertion_code}"
                raise ValidationError("PDB can not contain residues with insertion codes: " + res_id)


class JobSubmitForm(forms.Form):
    job_name = forms.CharField(label='Job name',
                               help_text='Enter your job name here',
                               max_length=100,
                               required=False)

    rec_pdb_file = forms.FileField(label='Protein *',
                                   help_text='Provide your docking target',
                                   required=True,
                                   validators=[validate_file, validate_pdb])

    rec_chain = forms.CharField(label='Chain ID *',
                                help_text='Enter name of the protein chain',
                                max_length=1,
                                required=True)

    lig_smiles = forms.CharField(label='Ligand SMILES *',
                                 help_text='Enter your ligand in SMILES format',
                                 required=True)

    pdb_exc = forms.CharField(label='PDB exclusions',
                              help_text='Enter a list of PDB codes, which will be excluded from homology search (comma separated)',
                              required=False)

    modeller_key = forms.CharField(label='MODELLER key',
                                   help_text='Enter your MODELLER key or obtain one from <a href=https://salilab.org/modeller/registration.html>here</a>',
                                   required=False)

    no_remodeling = forms.BooleanField(label='Do not remodel',
                                       help_text="Check if you don't want your protein to be remodeled based on chosen template",
                                       required=False)

    def clean_lig_smiles(self):  # TODO: add smiles validation
        smi = self.cleaned_data.get('lig_smiles')
        if len(smi.split()) != 1:
            self.add_error('lig_smiles', 'Only single ligand can be supplied at a time. Spaces and line-breaks are not allowed.')
        return smi

    def clean_pdb_exc(self):
        pdb_list = self.cleaned_data.get('pdb_exc')
        if not pdb_list:
            return []

        try:
            pdb_list = pdb_list.split(',')
            pdb_list = list(map(lambda x: x.strip().upper(), pdb_list))
        except Exception as e:
            self.add_error('pdb_exc', 'Wrong format of PDB exclusions')
            return []

        for pdb in pdb_list:
            if len(pdb) != 4:
                self.add_error('pdb_exc', 'Wrong PDB code: ' + pdb)
                return []

        return pdb_list

    def clean_rec_chain(self):
        if len(self.cleaned_data.get('rec_chain')) != 1:
            self.add_error('rec_chain', 'Chain name must be a single letter')

        return self.cleaned_data.get('rec_chain')

    def clean(self):
        cleaned_data = super().clean()
        if not cleaned_data.get('no_remodeling') and cleaned_data.get('modeller_key') == '':
            self.add_error('modeller_key', 'MODELLER key is required unless "Do not remodel" is checked')

        return cleaned_data


class SettingsForm(forms.Form):
    current_password = forms.CharField(label='Current password',
                                       max_length=100,
                                       required=True,
                                       help_text='Please enter your current password',
                                       widget=forms.PasswordInput(attrs={'class': 'form-control _password'}))

    password = forms.CharField(label='Password',
                               max_length=100,
                               required=True,
                               validators=[validate_password],
                               help_text=password_validators_help_text_html(),
                               widget=forms.PasswordInput(attrs={'class': 'form-control _password'}))


class AcademicEmailField(forms.EmailField):
    def validate(self, value):
        super().validate(value)
        # List based on Wikipedia (https://en.wikipedia.org/wiki/.edu_(second-level_domain) and https://en.wikipedia.org/wiki/.edu_(second-level_domain))
        tld_good = (r'\.edu$', r'\.edu\.[a-z]+$', r'\.ac\.[a-z]+$')
        if not any(re.match(r'^.*@.*'+tld, value) for tld in tld_good):
            raise ValidationError('Your e-mail must belong to an academic domain (.edu, .edu.*, .ac.*)')
        if User.objects.filter(email__exact=value).exists():
            raise ValidationError('Provided email is already in the database')


class SignUpForm(forms.Form):
    email = AcademicEmailField(label='E-mail *',
                               required=True,
                               max_length=1000,
                               help_text='Please, provide a valid academic e-mail address (.edu, .edu.*, .ac.*), we will use it to send you your password',
                               widget=forms.EmailInput(attrs={'class': 'form-control',
                                                              'style': 'width: 40ch'}))

    first_name = forms.CharField(label='First Name',
                                 max_length=1000,
                                 required=False,
                                 widget=forms.TextInput(attrs={'class': 'form-control',
                                                               'style': 'width: 40ch'}))

    last_name = forms.CharField(label='Last Name',
                                max_length=1000,
                                required=False,
                                widget=forms.TextInput(attrs={'class': 'form-control',
                                                              'style': 'width: 40ch'}))


def username_not_exists_validator(value):
    if not User.objects.filter(username__exact=value).exists():
        raise ValidationError('Sorry, we couldn\'t find the e-mail you provided')


class PasswordResetForm(forms.Form):
    username = forms.CharField(label='E-mail',
                               required=True,
                               max_length=1000,
                               help_text='Please, enter your e-mail address',
                               validators=[username_not_exists_validator],
                               widget=forms.TextInput(attrs={'class': 'form-control',
                                                             'style': 'width: 40ch'}))
