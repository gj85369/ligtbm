#define _POSIX_C_SOURCE 200809L
#include <getopt.h>
#include <jansson.h>
#include <stdbool.h>
#include <mol2/atom_group.h>
#include <mol2/gbsa.h>
#include <mol2/icharmm.h>
#include <mol2/json.h>
#include <mol2/minimize.h>
#include <mol2/pdb.h>
#include <mol2/prms.h>
#include <mol2/fitting.h>
#include <sampling/energy.h>
#include <sampling/transform.h>
#include <sampling/utils.h>

#define __FILENAME__ "main.c"


void __die(const char *message_format, va_list message_args) __attribute__((noreturn));
void __die(const char *message_format, va_list message_args)
{
	fprintf(stderr, "Error: ");
	vfprintf(stderr, message_format, message_args);
	fprintf(stderr, "\n");
	exit(EXIT_FAILURE);
}


void _die_if(bool is_error, const char *message_format, ...)
{
	va_list message_args;
	va_start (message_args, message_format);
	if (is_error) {
		__die(message_format, message_args);
	}
	va_end(message_args);
}


void _freeze_all_unfreeze_atom_list(struct mol_atom_group *ag, const struct mol_index_list atom_list)
{
	for (size_t i = 0; i < ag->natoms; i++) {
		ag->fixed[i] = true;
	}

	for (size_t i = 0; i < atom_list.size; i++) {
		const size_t j = atom_list.members[i];
		ag->fixed[j] = false;
	}
}


void _unfreeze_atoms_from_list_if_near(
	struct mol_atom_group *ag,
	const struct mol_index_list list_to_unfreeze,
	const double range,
	const struct mol_index_list list_reference)
{
	const double range_sq = range * range;
	for (size_t i = 0; i < list_to_unfreeze.size; i++) {
		const size_t atom_i = list_to_unfreeze.members[i];
		for (size_t j = 0; j < list_reference.size; j++) {
			const size_t atom_j = list_reference.members[j];
			const double dist_sq = MOL_VEC_EUCLIDEAN_DIST_SQ(ag->coords[atom_i], ag->coords[atom_j]);
			bool is_in_range = (dist_sq < range_sq);
			ag->fixed[atom_i] = !is_in_range;
		}
	}
}


struct detailed_energy {
	double genborn;
	double vdw_rep;
	double vdw_atr;
	double coulomb_shrt;
	double coulomb_long;
	double bonds;
	double angles;
	double impropers;
	double torsions;
	double restraints;
	double density;
};


struct restraints_params {
	struct mol_index_list atom_indices;
	struct mol_vector3 *coords;
	double k_spring;
};


struct density_params {
    struct mol_atom_group* ref_ag;
    bool no_grad;
    bool *mask;
    double weight;
    double radius;
};


struct potential_params {
	struct mol_atom_group *ag;
	struct agsetup *ags;
	struct acesetup *acs;
	struct mol_index_list lig_list;
	struct detailed_energy energy;
	struct restraints_params restraints;
	struct density_params density;
	bool enable_torsions;
	bool enable_electrostatics;
	bool enable_vdw;
	bool enable_density;
};


static void _zero_energies(struct detailed_energy *e)
{
	e->genborn = 0;
	e->vdw_rep = 0;
	e->vdw_atr = 0;
	e->coulomb_shrt = 0;
	e->coulomb_long = 0;
	e->bonds = 0;
	e->angles = 0;
	e->impropers = 0;
	e->torsions = 0;
	e->restraints = 0;
	e->density = 0;
}

static double _total_energy(const struct detailed_energy *e)
{
	return e->genborn + e->vdw_rep + e->vdw_atr + \
		e->coulomb_shrt + e->coulomb_long + \
		e->bonds + e->angles + e->impropers + e->torsions + \
		e->restraints + e->density;
}

static void _write_energies(
		const struct detailed_energy *e,
		const char *ofile)
{
	FILE *fout = fopen(ofile, "w");
	_die_if(fout == NULL,
		"Unable to open %s for writing", ofile);

	// It's a handmade JSON. I feel shame.
	fprintf(fout, "{");
#define _PRINT_TERM(name) do { fprintf(fout, "\"%s\": %.5f, ", #name, e->name); } while(0)
	_PRINT_TERM(genborn);
	_PRINT_TERM(vdw_rep);
	_PRINT_TERM(vdw_atr);
	_PRINT_TERM(coulomb_shrt);
	_PRINT_TERM(coulomb_long);
	_PRINT_TERM(bonds);
	_PRINT_TERM(angles);
	_PRINT_TERM(impropers);
	_PRINT_TERM(torsions);
	_PRINT_TERM(restraints);
	_PRINT_TERM(density);
#undef _PRINT_TERM
	double total = _total_energy(e);
	fprintf(fout, "\"total\": %.5f", total);
	fprintf(fout, "}\n");
	fclose(fout);
}


static void _set_beta_factors_for_output(
		struct mol_atom_group *ag,
                const struct mol_index_list list_rec,
                const struct mol_index_list list_lig,
                const struct mol_index_list list_restraints)
{
	// Mark beta: 0 - frozen atom; 1 - free receptor atom (at least on the last stage); 2 - free ligand atom; 3 - restrained ligand atom.
	for (size_t atom = 0; atom < ag->natoms; atom++)
		ag->B[atom] = 0;
	for (size_t i = 0; i < list_rec.size; i++) {
		size_t atom = list_rec.members[i];
		if (!ag->fixed[atom])
			ag->B[atom] = 1;
	}
	for (size_t i = 0; i < list_lig.size; i++) {
		size_t atom = list_lig.members[i];
		ag->B[atom] = 2;
	}
	for (size_t i = 0; i < list_restraints.size; i++) {
		size_t atom = list_restraints.members[i];
		ag->B[atom] = 3;
	}
}


// The \p atom_list contains the atoms in the atom group whose coordinates will be copied into \p array.
static void _coords_to_array_atom_list(
	struct mol_vector3 *dst,
	const struct mol_atom_group *ag,
	const struct mol_index_list atom_list)
{
	for (size_t j = 0; j < atom_list.size; j++) {
		size_t i = atom_list.members[j];
		MOL_VEC_COPY(dst[j], ag->coords[i]);
	}
}


static void _find_common_atoms_by_name(
		const struct mol_atom_group *ag_1,
		const struct mol_atom_group *ag_2,
		const struct mol_index_list list_to_use_in_ag_2,
		struct mol_index_list *list_common_1,
		struct mol_index_list *list_common_2)
{
#define MAX_ATOM_NAME_LEN 20  // Chosen semi-arbitrary
	size_t count = 0;
	size_t max_atoms = MAX(ag_1->natoms, ag_2->natoms);
	list_common_1->members = calloc(max_atoms, sizeof(size_t));
	list_common_2->members = calloc(max_atoms, sizeof(size_t));
	char name_i[MAX_ATOM_NAME_LEN], name_j[MAX_ATOM_NAME_LEN];

	for (size_t atom_i = 0; atom_i < ag_1->natoms; atom_i++) {
		smp_str_trim_and_copy(name_i, ag_1->atom_name[atom_i], MAX_ATOM_NAME_LEN);
		for (size_t j = 0; j < list_to_use_in_ag_2.size; j++) {
			size_t atom_j = list_to_use_in_ag_2.members[j];
			smp_str_trim_and_copy(name_j, ag_2->atom_name[atom_j], MAX_ATOM_NAME_LEN);
			if (strcmp(name_i, name_j) == 0) {
				list_common_1->members[count] = atom_i;
				list_common_2->members[count] = atom_j;
				count++;
				break;
			}
		}
	}

	list_common_1->members = realloc(list_common_1->members, count * sizeof(size_t));
	list_common_2->members = realloc(list_common_2->members, count * sizeof(size_t));
	list_common_1->size = count;
	list_common_2->size = count;
#undef MAX_ATOM_NAME_LEN
}

double potential(struct potential_params *pot_prms)
{
	static const double weight_vdw_atr = 1.0;
	static const double weight_vdw_rep = 1.0;
	static const double weight_coulomb_short = 1.0;
	static const double weight_coulomb_long = 1.0;

	struct detailed_energy *energies = &(pot_prms->energy);

	bool nblist_updated = check_clusterupdate(pot_prms->ag, pot_prms->ags);
	if (nblist_updated) {
		ace_updatenblst(pot_prms->ags, pot_prms->acs);
	}

	_zero_energies(energies);
	mol_zero_gradients(pot_prms->ag);

	if (pot_prms->enable_vdw) {
		smp_energy_vdw_rep_line(pot_prms->ag, &(energies->vdw_rep), &(energies->vdw_atr), pot_prms->ags->nblst,
		                        weight_vdw_rep, weight_vdw_atr);
		smp_energy_vdw03_rep_line(pot_prms->ag, &(energies->vdw_rep), &(energies->vdw_atr),
		                          pot_prms->ags->nblst, pot_prms->ags->nf03, pot_prms->ags->listf03,
		                          weight_vdw_rep, weight_vdw_atr);
	}
							  
	if (pot_prms->enable_electrostatics) {
		smp_energy_coul_rdie_shift(pot_prms->ag, &(energies->coulomb_shrt), &(energies->coulomb_long),
		                           pot_prms->ags->nblst,
		                           weight_coulomb_short, weight_coulomb_long);
		aceeng(pot_prms->ag, &(energies->genborn), pot_prms->acs, pot_prms->ags);
	}

	if (pot_prms->enable_density) {
	    const struct mol_fitting_params den_prms = {.radius = pot_prms->density.radius,
	                                                .mask = pot_prms->density.mask,
	                                                .normalize = true};

	    struct mol_vector3* grads = pot_prms->ag->gradients;
	    if (pot_prms->density.no_grad) {
            pot_prms->ag->gradients = NULL;
	    }
        energies->density = mol_fitting_score(pot_prms->ag,
                                              pot_prms->density.ref_ag,
                                              &den_prms,
                                              pot_prms->density.weight);
        pot_prms->ag->gradients = grads;
	}

	smp_energy_restraints_fixedpoint(pot_prms->ag, &(energies->restraints),
		pot_prms->restraints.atom_indices,
		pot_prms->restraints.coords,
		pot_prms->restraints.k_spring);

	beng(pot_prms->ag, &(energies->bonds));
	aeng(pot_prms->ag, &(energies->angles));
	ieng(pot_prms->ag, &(energies->impropers));
	if (pot_prms->enable_torsions)
		teng(pot_prms->ag, &(energies->torsions));

    if (pot_prms->density.no_grad) {
        return _total_energy(energies) - energies->density;
    }
    return _total_energy(energies);
}

// Passed as a function pointer to the LBFGS
lbfgsfloatval_t potential_wrapper(
		void *prms,
		const double *inp,
		double *grad,
		const int n,
		__attribute__((unused)) const double step)
{
	struct potential_params *pot_prms = (struct potential_params*) prms;
	struct mol_atom_group *ag = pot_prms->ag;

	if (inp != NULL) {
		_die_if(n != ((int) ag->active_atoms->size) * 3,
			"Mismatch in in vector length");
		mol_atom_group_set_actives(ag, inp);
	}

	double energy = potential(pot_prms);

	if (grad != NULL) {
		for (int i = 0; i < n / 3; i++) {
			struct mol_vector3 neg_grad;
			MOL_VEC_MULT_SCALAR(neg_grad, ag->gradients[ag->active_atoms->members[i]], -1);
			SMP_VEC_TO_ARRAY(grad, i, neg_grad);
		}
	}
	return energy;
}


static void run_all_atom_min_simple(
		struct potential_params *pot_prms,
		const double cutoff,
		const double restraint_spring,
		const struct mol_index_list receptor_hydrogens)
{
	static const double min_tol = 1e-6;
	struct mol_atom_group *ag = pot_prms->ag;
	struct agsetup *ags = pot_prms->ags;
	struct acesetup *acs = pot_prms->acs;
	struct mol_index_list lig_list = pot_prms->lig_list;

	// Unfreeze receptor hydrogens near ligand
	_freeze_all_unfreeze_atom_list(ag, lig_list);
	_unfreeze_atoms_from_list_if_near(ag, receptor_hydrogens, cutoff, lig_list);

	mol_fixed_update_active_lists(ag);
	update_nblst(ag, ags);
	ace_fixedupdate(ag, ags, acs);
	ace_updatenblst(ags, acs);

	pot_prms->enable_electrostatics = true;
	pot_prms->enable_vdw = true;
	pot_prms->enable_torsions = true;
	pot_prms->enable_density = false;

	pot_prms->restraints.k_spring = 0;
	mol_minimize_ag(MOL_LBFGS, 500, min_tol, ag, (void *) pot_prms, potential_wrapper);

	pot_prms->restraints.k_spring = restraint_spring;
	mol_minimize_ag(MOL_LBFGS, 500, min_tol, ag, (void *) pot_prms, potential_wrapper);

	pot_prms->restraints.k_spring = 0;
	pot_prms->enable_density = true;
	mol_minimize_ag(MOL_LBFGS, 500, min_tol, ag, (void *) pot_prms, potential_wrapper);

	pot_prms->enable_density = false;
	mol_minimize_ag(MOL_LBFGS, 500, min_tol, ag, (void *) pot_prms, potential_wrapper);
}



static void run_all_atom_min_multistage(
		struct potential_params *pot_prms,
		const double cutoff,
		const double restraint_spring,
		const struct mol_index_list receptor_hydrogens,
		const struct mol_index_list receptor_all)
{
	static const double min_tol = 1E-6;
	struct mol_atom_group *ag = pot_prms->ag;
	struct agsetup *ags = pot_prms->ags;
	struct acesetup *acs = pot_prms->acs;
	struct mol_index_list lig_list = pot_prms->lig_list;

	// Unfreeze receptor hydrogens near ligand
	_freeze_all_unfreeze_atom_list(ag, lig_list);
	_unfreeze_atoms_from_list_if_near(ag, receptor_hydrogens, cutoff, lig_list);
	mol_fixed_update_active_lists(ag);
	update_nblst(ag, ags);
	ace_fixedupdate(ag, ags, acs);
	ace_updatenblst(ags, acs);

	pot_prms->enable_torsions = false;
	pot_prms->enable_electrostatics = false;
	pot_prms->enable_vdw = false;
	pot_prms->enable_density = false;
	pot_prms->restraints.k_spring = restraint_spring;
	mol_minimize_ag(MOL_LBFGS, 500, min_tol, ag, (void *) pot_prms, potential_wrapper);

	pot_prms->enable_torsions = true;
	mol_minimize_ag(MOL_LBFGS, 500, min_tol, ag, (void *) pot_prms, potential_wrapper);

	_unfreeze_atoms_from_list_if_near(ag, receptor_all, cutoff, lig_list);
	mol_fixed_update_active_lists(ag);
	update_nblst(ag, ags);
	ace_fixedupdate(ag, ags, acs);
	ace_updatenblst(ags, acs);

	pot_prms->enable_electrostatics = true;
	pot_prms->enable_vdw = true;
	mol_minimize_ag(MOL_LBFGS, 500, min_tol, ag, (void *) pot_prms, potential_wrapper);

	pot_prms->restraints.k_spring = 0;
	mol_minimize_ag(MOL_LBFGS, 500, min_tol, ag, (void *) pot_prms, potential_wrapper);
}


void print_usage(const char *self_name)
{
	fprintf(stderr, "RMIN version %s\n\n", _GIT_VERSION_);
	fprintf(stderr, "USAGE: %s " \
		 "<rec.pdb> <rec.psf> " \
		 "<lig.json> <lig_ref.pdb> " \
		 "<forcefield.prm> <topology.rtf> <libmol.params> " \
		 "<out_min.pdb> <out_min.dat> [options]\n", self_name);
	fprintf(stderr, "   rec.pdb           - Receptor PDB file.\n");
	fprintf(stderr, "   rec.psf           - Receptor PSF file.\n");
	fprintf(stderr, "   lig.json          - Ligand structure in libmol JSON format.\n");
	fprintf(stderr, "   lig_ref.pdb       - Ligand reference PDB file (for RMSD and geometric restraints).\n");
	fprintf(stderr, "   forcefield.prm    - CHARMM-style forcefield parameters.\n");
	fprintf(stderr, "   topology.rtf      - CHARMM-style topology file.\n");
	fprintf(stderr, "   libmol.params     - libmol-style atom type parameters.\n");
	fprintf(stderr, "   out_min.pdb       - PDB file to which lowest-energy structure should be written.\n");
	fprintf(stderr, "   out_min.dat       - File to which lowest-energy values should be written.\n");
	fprintf(stderr, " Options:\n");
	fprintf(stderr, "  -s  --scale=<n>    - Scale default restraints spring constant. Default: 1.\n");
	fprintf(stderr, "  -d  --density=<n>  - Scale default density constant. Default: 1.\n");
	fprintf(stderr, "  -c  --cutoff=<n>   - Receptor cutoff. Distance from the ligand center for which the receptor atoms are unfrozen. Default: 5 A.\n");
	fprintf(stderr, "  -p  --protocol=<n> - Minimization protocol. 1 - Simple, 2 - Multistage. Default: 1\n");
	fprintf(stderr, "  -h  --help         - Print this message and quit.\n");
}


int main(int argc, char **argv)
{
	if (argc < 10) {
		print_usage(argv[0]);
		exit(EXIT_FAILURE);
	}
	// Default values
	double cutoff_around_ligand = 5.0; // Angstrom
	double restraints_k_spring = 1.0;
	double density_weight = 1.0;
	double density_radius = 2.0;
	int protocol = 1;

	// Buffers
	const char *rec_pdb_file = strdup(argv[1]);
	const char *rec_psf_file = strdup(argv[2]);
	const char *lig_json_file = strdup(argv[3]);
	const char *pdb_file_restraints = strdup(argv[4]);
	const char *ff_prm_file = strdup(argv[5]);
	const char *ff_rtf_file = strdup(argv[6]);
	const char *ff_libmol_file = strdup(argv[7]);
	const char *out_pdb_file = strdup(argv[8]);
	const char *out_dat_file = strdup(argv[9]);

	static const struct option long_options[] = {
		{"spring", required_argument, NULL, 's'},
		{"density", required_argument, NULL, 'd'},
		{"cutoff", required_argument, NULL, 'c'},
		{"protocol", required_argument, NULL, 'p'},
		{"help", no_argument, NULL, 'h'},
		{0, 0, 0, 0}
	};

	while(1) {
		int option_index;
		int c = getopt_long(argc, argv, "d:s:c:p:h", long_options, &option_index);
		if (c == -1)
			break;
		switch(c) {
			case 'h':
				print_usage(argv[0]);
				exit(EXIT_SUCCESS);
			case 's':
				restraints_k_spring = atol(optarg);
				break;
			case 'd':
				density_weight = atof(optarg);
				break;
			case 'c':
				cutoff_around_ligand = atof(optarg);
				break;
			case 'p':
				protocol = atoi(optarg);
				break;
		}
	}

	printf("Core parameters:\n");
	printf(" Receptor PDB: %s\n          PSF: %s\n", rec_pdb_file, rec_psf_file);
	printf(" Ligand JSON: %s\n        PDB restraints: %s\n", lig_json_file, pdb_file_restraints);
	printf(" Forcefield parameters: %s\n            geometry: %s\n            atom types: %s\n", ff_prm_file, ff_rtf_file, ff_libmol_file);
	printf(" Output final structure PDB: %s\n        final energy DAT: %s\n", out_pdb_file, out_dat_file);
	printf(" Cutoff around ligand: %f A\n Restraints spring: %f\n", cutoff_around_ligand, restraints_k_spring);
	printf(" Density weight: %f\n Protocol number: %d\n", density_weight, protocol);

	struct mol_prms *atomprm = mol_prms_read(ff_libmol_file);
	_die_if(atomprm == NULL, "Unable to load libmol parameters");

	// Load receptor
	struct mol_atom_group *ag_rec = mol_read_pdb(rec_pdb_file);
	_die_if(ag_rec == NULL, "Unable to load the receptor");
	bool read_ok = mol_atom_group_read_geometry(ag_rec, rec_psf_file, ff_prm_file, ff_rtf_file);
	_die_if(!read_ok, "Unable to load geometry (psf and forcefield) for the receptor");
	mol_atom_group_add_prms(ag_rec, atomprm);

	// Set up list of receptor atoms in complex that will be created later
	struct mol_index_list list_rec;
	list_rec.size = ag_rec->natoms;
	list_rec.members = calloc(list_rec.size, sizeof(size_t));
	for (size_t i = 0; i < ag_rec->natoms; i++)
		list_rec.members[i] = i;

	// Load ligand
	struct mol_atom_group *ag_lig = mol_read_json(lig_json_file);
	_die_if(ag_lig == NULL, "Unable to load the ligand");

	// Set up list of ligand atoms in complex that will be created later
	struct mol_index_list list_lig;
	list_lig.size = ag_lig->natoms;
	list_lig.members = calloc(list_lig.size, sizeof(size_t));
	for (size_t i = 0; i < list_lig.size; i++)
		list_lig.members[i] = ag_rec->natoms + i;

	// Load restraints
	struct mol_atom_group *ag_reference = mol_read_pdb(pdb_file_restraints);
	_die_if(ag_reference == NULL, "Unable to load the restraint file");

	// Create complex
	struct mol_atom_group *ag = mol_atom_group_join(ag_rec, ag_lig);
	mol_atom_group_create_residue_list(ag);
	mol_fixed_init(ag);

	// Setup non-bonded list
	struct agsetup ags;
	init_nblst(ag, &ags);
	update_nblst(ag, &ags);

	// Setup ACE solvent model
	struct acesetup acs;
	acs.efac = 0.5; // Use the same efac as in David's minimization code
	ace_ini(ag, &acs);
	ace_fixedupdate(ag, &ags, &acs);
	ace_updatenblst(&ags, &acs);

	// For all-atom minimization find hydrogen receptors
	struct mol_index_list list_rec_hydrogens;
	smp_get_atom_list_by_element(&list_rec_hydrogens, ag_rec, "H", true);

	// Find atoms in ligand for which reference constraints
	struct mol_index_list list_restrained_in_ag;
	struct mol_index_list list_restrained_in_ref;
	_find_common_atoms_by_name(ag_reference, ag, list_lig, &list_restrained_in_ref, &list_restrained_in_ag);
	_die_if(list_restrained_in_ag.size == 0,
		"No common atoms between ligand and restraints references!");
	struct mol_vector3 *restraints_coords = calloc(list_restrained_in_ref.size, sizeof(struct mol_vector3));
	_coords_to_array_atom_list(restraints_coords, ag_reference, list_restrained_in_ref);

	// Initialize structure with energy calculation parameters
	// Set up pointers to current ag, ags, and all other needed data.
	struct potential_params pot_prms;
	pot_prms.ag = ag;
	pot_prms.ags = &ags;
	pot_prms.acs = &acs;
	pot_prms.lig_list = list_lig;
	pot_prms.restraints.k_spring = restraints_k_spring;
	pot_prms.restraints.atom_indices = list_restrained_in_ag;
	pot_prms.restraints.coords = restraints_coords;
	pot_prms.density.ref_ag = ag_reference;
	pot_prms.density.weight = density_weight;
	pot_prms.density.radius = density_radius;
	_zero_energies(&pot_prms.energy);

	// Setup density params
	bool *density_mask = calloc(ag->natoms, sizeof(bool));
	for (size_t i = 0; i < ag_rec->natoms; i++) {
		density_mask[i] = false;
	}
	for (size_t i = ag_rec->natoms; i < ag->natoms; i++) {
		density_mask[i] = true;
	}
	pot_prms.density.mask = density_mask;
	pot_prms.density.no_grad = false;

	// All-atom minimization
	switch (protocol) {
		case 1:
			run_all_atom_min_simple((void *) (&pot_prms), cutoff_around_ligand, restraints_k_spring, list_rec_hydrogens);
			break;
		case 2:
			run_all_atom_min_multistage((void *) (&pot_prms), cutoff_around_ligand, restraints_k_spring, list_rec_hydrogens, list_rec);
			break;
		default:
			fprintf(stderr, "[error] Wrong protocol number\n");
			exit(EXIT_FAILURE);
	}

	// Save low energy structure
	_set_beta_factors_for_output(ag, list_rec, list_lig, list_restrained_in_ag);
	mol_write_pdb(out_pdb_file, ag);
	_write_energies(&(pot_prms.energy), out_dat_file);

	free(density_mask);

	return 0;
}
