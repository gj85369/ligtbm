function fill_stage(stage, data, model_id) {
    let blob = new Blob([data], {type: 'text/plain'});
    stage.removeAllComponents();
    stage.loadFile(blob, {ext: 'pdb', name: 'mol'}).then(function (c) {
        c.addRepresentation("cartoon", {
            sele: "not :X",
            scale: 1.0,
            smoothSheet: true
        });
        //c.addRepresentation("line", {sele: "not :X"});
        c.addRepresentation("licorice", {sele: ":X"});
        c.autoView(":X");
        c.setSelection("/".concat(String(model_id)));
    });
}

$(document).ready(function () {
    let slider = document.getElementById("myRange");
    let output = document.getElementById("demo");
    output.innerHTML = "Model " + slider.value;

    let stage = new NGL.Stage("viewport", {backgroundColor: "white"});
    let models;
    $.get($("#models_link").attr("href"), function (data) {
        models = data;
        //console.log(data);
    }).then(function () {
        let last_model_id = Math.min(models.split('ENDMDL').length - 1, 10);
        //console.log(last_model_id);
        $("#myRange").attr("max", last_model_id);

        // Fill the stage with the first model
        fill_stage(stage, models, 0);
    });

    // Update slider value each time the handle is moved
    slider.oninput = function () {
        output.innerHTML = "Model " + slider.value;
        let model_id = parseInt(this.value) - 1;
        //fill_stage(stage, models, model_id);
        let comp = stage.getComponentsByName('mol');
        //console.log(comp.list[0]);
        comp.list[0].setSelection("/".concat(String(model_id)));
    }
});