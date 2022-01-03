console.log("active_site_for_chain");
console.log(active_site_for_chain);
console.log("all_ligand_binding_for_chain");
console.log(all_ligand_binding_for_chain);
console.log("interfacing_for_chain");
console.log(interfacing_for_chain);
console.log("dna_binding_for_chain");
console.log(dna_binding_for_chain);
console.log("intrachain_disulphide_list");
console.log(intrachain_disulphide_list);
salt_bridges = Array.from(new Set(salt_bridges));
console.log("salt_bridges");
console.log(salt_bridges);
console.log("");
console.log("");
console.log("");
console.log("prediction_per_pos");
console.log(prediction_per_pos);
console.log("prediction_per_mut");
console.log(prediction_per_mut);
console.log("coevolution_pos");
console.log(coevolution_pos);
console.log("additional_mappings");
console.log(additional_mappings);

var res_keys = Object.keys(prediction_per_pos);
var sorted_res_keys = res_keys.sort(function(a, b) {
    var anum = parseInt(a.substr(1, a.length-1));
    var bnum = parseInt(b.substr(1, b.length-1));
    return anum - bnum;
});
var total_res_num = Object.keys(prediction_per_pos).length;
var prediction_per_mut_f = Object.values(prediction_per_mut).filter(function(a){
    return a[1] !== "None";
});
var total_mut_num = Object.keys(prediction_per_mut_f).length;
var box_per_res = 35;
var total_width_general = total_res_num * box_per_res;
var width_80_perc = window.innerWidth * .8;
var width_10_perc = window.innerWidth * .1;
var width_65_perc = window.innerWidth * .65;
var total_width_general_s = width_80_perc/2;

var line_plot_x = 30;
var line_plot_sx = 30;
var line_plot_y = 15;

var mutation_list = [
    "A","C","D","E","F",
    "G","H","I","K","L",
    "M","N","P","Q","R",
    "S","T","V","W","Y"
];

var threeToOne = {
    "ALA": "A", "MET": "M",
    "CYS": "C", "ASN": "N",
    "ASP": "D", "PRO": "P",
    "GLU": "E", "GLN": "Q",
    "PHE": "F", "ARG": "R",
    "GLY": "G", "SER": "S",
    "HIS": "H", "THR": "T",
    "ILE": "I", "VAL": "V",
    "LYS": "K", "TRP": "W",
    "LEU": "L", "TYR": "Y"
}
var oneToThree = {
    "A":"ALA", "M":"MET",
    "C":"CYS", "N":"ASN",
    "D":"ASP", "P":"PRO",
    "E":"GLU", "Q":"GLN",
    "F":"PHE", "R":"ARG",
    "G":"GLY", "S":"SER",
    "H":"HIS", "T":"THR",
    "I":"ILE", "V":"VAL",
    "K":"LYS", "W":"TRP",
    "L":"LEU", "Y":"TYR"
}


function formatVal(num, after_comma) {
    let ten_basis = Math.pow(10, after_comma)
    return Math.round( num * ten_basis + Number.EPSILON ) / ten_basis;
}

var prediction_means = []
var prediction_amplitudes = [];
var all_prediction_values = [];
var yScale;
function generate_dsm_lineplot(svg_lineplot, legend_svg) {
    // get means and amplitude of mutation effects
    // prediction_per_pos
    prediction_means = [];
    prediction_amplitudes = [];
    all_prediction_values = [];
    let i = 0.5;
    for (const reskey of sorted_res_keys) {
        var values = prediction_per_pos[reskey];
        var prediction_values = values.map(function(a){
            return a[1];
        })
        prediction_values = prediction_values.filter(function(a) {
            return a !== "None";
        });
        prediction_means.push({"pos_id": reskey,"x": (i*box_per_res)+line_plot_x,"raw_y":d3.mean(prediction_values)});
        prediction_amplitudes.push({"x": (i*box_per_res)+line_plot_x,"raw_y_array":d3.extent(prediction_values)});
        all_prediction_values.push(...prediction_values);
        i += 1;
    }
    let y_extent = d3.extent(all_prediction_values);

    yScale = d3.scaleLinear()
    .domain(y_extent)
    .range([200-line_plot_y, line_plot_y]);

    // legend_svg.append("g")
    // .attr("transform", "translate(0," + height + ")")
    // .call(d3.axisBottom(xScale));
   
    // svg_lineplot_yaxes.append("g")
    // legend_svg.append("g")
    // .attr("transform", "translate("+(width_10_perc-5)+",-15)")
    // .call(d3.axisLeft(yScale));

    for (let i = 0; i < prediction_means.length; i++) {
        prediction_means[i]["y"] = yScale(prediction_means[i]["raw_y"]);
        prediction_amplitudes[i]["y_array"] = [yScale(prediction_amplitudes[i]["raw_y_array"][0]),yScale(prediction_amplitudes[i]["raw_y_array"][1])];
    }

    var line = d3.line()
        .x(function(d) { return d.x; }) 
        .y(function(d) { return d.y; }) 
        .curve(d3.curveMonotoneX);

    var wt_line_datum = [
        {"x": 0.0, "y":yScale(0.0)},
        {"x": total_width_general+(line_plot_x*2), "y":yScale(0.0)},
    ]
    svg_lineplot.append("path")
        .datum(wt_line_datum) 
        .attr("class", "line") 
        // .attr("transform", "translate(" + 100 + "," + 100 + ")")
        .attr("d", line)
        .style("fill", "none")
        .style("stroke", "grey")
        .style("stroke-width", "2")
        .on("mouseover", function(d) {		
            d3.select(this).style("stroke", "red");
            d3.select(this).style("stroke-width", "4");
            d3.select("#pagetooltip")
            .transition()
            .duration(200)
            .style("opacity", .9);
            d3.select("#pagetooltip")
            .html(`Wild type line`)
            .style("left", (d3.event.pageX) + "px")		
            .style("top", (d3.event.pageY - 28) + "px");
        })					
        .on("mouseout", function(d) {	
            d3.select(this).style("stroke", "grey");
            d3.select(this).style("stroke-width", "2");
            d3.select("#pagetooltip")
            .transition()		
            .duration(500)		
            .style("opacity", 0);
        });
    
    //FULL LINE PLOT
    svg_lineplot.append("path")
        .datum(prediction_means) 
        .attr("class", "line") 
        // .attr("transform", "translate(" + 100 + "," + 100 + ")")
        .attr("d", line)
        .style("fill", "none")
        .style("stroke", "black")
        .style("stroke-width", "2");
    
    //ERROR BAR PLOTTING
    svg_lineplot.append("g").selectAll(".ampl_bars")
    .data(prediction_amplitudes)
    .enter()
    .append("path")
    .attr("class", "ampl_bars")
    .attr("d", function(d){
        return line([{"x":d.x, "y":d.y_array[0]},{"x":d.x, "y":d.y_array[1]}])
    })
    .style("fill", "none")
    .style("stroke", "black")
    .style("stroke-width", "2");
    
    //MEAN VAL CIRCLE
    svg_lineplot.append("g").selectAll(".mean_circle")
    .data(prediction_means)
    .enter()
    .append("circle")
    .attr("id", function(d){ return d.pos_id})
    .attr("class", "mean_circle")
    .attr("cx", function(d){ return d.x})
    .attr("cy", function(d){ return d.y})
    .attr("r", 3)
    .attr("fill", "black")
    .on("mouseover", function(d) {		
        d3.select(this).attr("r", 6);
        d3.select(this).attr("fill", "red");
        d3.select("#pagetooltip")
        .transition()
        .duration(200)
        .style("opacity", .9);
        d3.select("#pagetooltip")
        .html(`Pos: ${d.pos_id}  Mean Score: ${formatVal(d.raw_y, 3)}`)
        .style("left", (d3.event.pageX) + "px")		
        .style("top", (d3.event.pageY - 28) + "px");
    })					
    .on("mouseout", function(d) {	
        d3.select(this).attr("r", 3);
        d3.select(this).attr("fill", "black");
        d3.select("#pagetooltip")
        .transition()		
        .duration(500)		
        .style("opacity", 0);
    })
    .on("click", function(d) {
        //TODO
        // add y axis values
        // add line for wild type
        // 
        // add sequence on x axis
        generate_position_view(d.pos_id);
    });

    // apply yScale to prediction_means and prediction_amplitudes
    // svg_lineplot
}

var rectangles_array = [];
var dsm_heatmap_height = 400;
function generate_dsm_heatmap(svg_rectplot, legend_svg) {
    // svg_rectplot
    rectangles_array = [];
    let i = 0.5;

    var yScale2 = d3.scaleLinear()
    .domain([mutation_list.length-1, 0])
    .range([dsm_heatmap_height-line_plot_y, line_plot_y]);

    var all_prediction_values2 = [];
    for (const reskey of sorted_res_keys) {
        var values_pos = prediction_per_pos[reskey];
        var prediction_values_pos = values_pos.map(function(a){
            return a[1];
        })
        prediction_values_pos = prediction_values_pos.filter(function(a) {
            return a !== "None";
        });
        all_prediction_values2.push(...prediction_values_pos);
        for (const mut of mutation_list) {
            var mut_key = reskey+mut;
            var values = prediction_per_mut[mut_key];
            rectangles_array.push({
                "respos":reskey,"mutpos":mut_key,
                "x": (i*box_per_res)+(line_plot_x/2), "y":yScale2(mutation_list.indexOf(mut)),
                "shap":values[2], "c":values[1],
            });
        }
        i += 1;
    }
    let y_extent2 = d3.extent(all_prediction_values2);
    var cScale = d3.scaleLinear()
    .domain([y_extent2[0], 0.0, y_extent2[1]])
    .range(["#4d52ff", "white", "#ff4d50"]);
    rectangles_array = rectangles_array.map(function(a){
        a.cl = cScale(a.c);
        return a;
    });
    console.log("rectangles_array");
    console.log(rectangles_array);

    svg_rectplot.selectAll(".heatmap_rect")
    .data(rectangles_array)
    .enter()
    .append("rect")
    .attr("class", "heatmap_rect")
    .attr("x", function(d){return d.x})
    .attr("y", function(d){return d.y})
    .attr("width", line_plot_x+4)
    .attr("height", yScale2(1)-yScale2(0))
    .attr("fill", function(d){
        if (d.c !== "None") {
            return d.cl;
        }
        return "grey";
    })
    .attr("stroke", "black")
    .attr("stroke-width", "")
    .on("mouseover", function(d) {		
        d3.select(this).attr("stroke-width", 4);
        d3.select("#pagetooltip")
        .transition()
        .duration(200)
        .style("opacity", .9);
        d3.select("#pagetooltip")
        .html(`Pos: ${d.respos}  Mut:${d.mutpos}  Score: ${formatVal(d.c, 3)}`)
        .style("left", (d3.event.pageX) + "px")		
        .style("top", (d3.event.pageY - 28) + "px");
        if (d.c !== "None") {
            d3.select("#svg_lineplot")
            .append("circle")
            .attr("id", "show_mut_circle")
            .attr("cx", function(){ return d.x+(line_plot_x/2) })
            .attr("cy", function(){ return yScale(d.c)})
            .attr("r", 4)
            .attr("fill", "red");
        }
    })					
    .on("mouseout", function(d) {	
        d3.select(this).attr("stroke-width", "");
        d3.select("#pagetooltip")
        .transition()		
        .duration(500)		
        .style("opacity", 0);
        if (d.c !== "None") {
            d3.select("#show_mut_circle").remove();
        }
    })
    .on("click", function(d) {
        //TODO
        // add y axis values
        // add line for wild type
        // 
        // add sequence on x axis
        generate_position_view(d.respos);
        generate_mutation_view(d.mutpos);
    });

    var legendCol = d3.legendColor()
    .scale(cScale)
    .cells(8);

    legend_svg.append("g")
    // .attr("transform", "translate(500,10)")
    .call(legendCol);


    legend_svg.selectAll(".mut_aas")
    .data(mutation_list)
    .enter()
    .append("text")
    .attr("class", "mut_aas")
    .attr("x", width_10_perc-15)
    .attr("y", function(d, i) {
        return yScale2(i)+2;
    })
    .text(function(d) {return d});

    // let legend_color_scale = d3.scaleBand()
    // .domain([0, mutation_list.length-1])
    // .range([(dsm_heatmap_height+line_plot_y), line_plot_y]);
}

function generate_dsm_view() {
    let d3_dsm_view = d3.select("#dsm_view");
    d3_dsm_view.html("");
    d3_dsm_view.style("text-align", "center")
    d3_dsm_view.style("max-height", width_80_perc+"px")
    d3_dsm_view.style("margin-left", width_10_perc+"px")
    d3_dsm_view.style("margin-right", width_10_perc+"px")

    d3_dsm_view.append("h1")
    .html("Main View")
    ;
    d3_dsm_view.append("h2")
    .html(`UniProt Code: ${uniprot_code} - PDB Chain: ${pdb_chain.replace(":", " ")}`)
    ;
    d3_dsm_view.append("h2")
    .html(`N. of Positions: ${total_res_num}`)
    ;
    d3_dsm_view.append("h2")
    .html(`N. of Mutations: ${total_mut_num}`)
    ;

    var svg_axes_containers = d3_dsm_view
    .append("div");

    var div_for_axes = svg_axes_containers.append("div")
    .attr("height", ((200+line_plot_y*2)+(dsm_heatmap_height+line_plot_y*2))+"px")
    .attr("width", width_10_perc+"px")
    .style("max-width", width_10_perc+"px")
    .style("display", "inline-block");

    var legend1_svg = div_for_axes.append("svg")
    .attr("height", ((200+line_plot_y*2)-8)+"px")
    .attr("width", width_10_perc+"px");

    var legend2_svg = div_for_axes.append("svg")
    .attr("height", ((dsm_heatmap_height+line_plot_y*2)+8)+"px")
    .attr("width", width_10_perc+"px");


    var svg_containers = svg_axes_containers
    .append("div")
    .attr("width", width_65_perc+"px")
    .style("max-width", width_65_perc+"px")
    .style("display", "inline-block")
    .style("overflow-x", "scroll");

    var svg_lineplot = svg_containers
    .append("svg")
    .attr("id", "svg_lineplot")
    .attr("width", total_width_general+line_plot_x*2)
    // .attr("width", total_width_general)
    .attr("height", (200+line_plot_y*2)+"px")
    ;
    generate_dsm_lineplot(svg_lineplot, legend1_svg);

    var svg_rectplot = svg_containers
    .append("svg")
    .attr("id", "svg_rectplot")
    .attr("width", total_width_general+line_plot_x*2)
    .attr("height", (dsm_heatmap_height+line_plot_y*2)+"px")
    ;
    generate_dsm_heatmap(svg_rectplot, legend2_svg);
}
generate_dsm_view();

function generate_position_view(respos) {
    let d3_position_view = d3.select("#position_view");
    let content = d3_position_view.html();
    if (!content) {
        d3_position_view.style("text-align", "center");
        d3_position_view.style("max-height", width_80_perc+"px");
        d3_position_view.style("margin-left", width_10_perc+"px");
        d3_position_view.style("margin-right", width_10_perc+"px");
        d3_position_view.append("br");
        d3_position_view.append("br");
        d3_position_view.append("hr");
        d3_position_view.append("h1")
        .html("Position Info")
        ;
        d3_position_view.append("div")
        .attr("id", "pos_info_div");
    } else {
        console.log("exists");
        d3.select("#pos_info_div").html("");
    }
    let resp_num = respos.substr(1, respos.length);
    let trespos = oneToThree[respos[0]] + resp_num;
    let unipresp_num = uniprot_residue_mapping[trespos];
    var trespos_unip  = oneToThree[respos[0]] + unipresp_num;
    let prediction_means_sorted = prediction_means.sort(function(a, b) {
        return b.raw_y - a.raw_y;
    })
    let prediction_means_ids = prediction_means_sorted.map(function(a) {return a.pos_id});
    let pos_rank = prediction_means_ids.indexOf(respos);
    let mean_score_val = prediction_means_sorted[pos_rank].raw_y;
    d3.select("#pos_info_div").append("h2")
    .html(`Position (PDB): ${trespos}`)
    ;

    d3.select("#pos_info_div").append("h2")
    .html(`Position(UniProt): ${trespos_unip}`)
    ;
    
    d3.select("#pos_info_div").append("h2")
    .html(`Mean Score: ${formatVal(mean_score_val, 3)}`)
    ;

    d3.select("#pos_info_div").append("h2")
    .html(`Mean Score Rank: ${pos_rank+1} of ${prediction_means_ids.length}`)
    .append("span")
    .html(' (Higher → More deleterious)').style("font-size", "16px");
    ;
    
    d3.select("#pos_info_div").append("br");

    // prediction_means
    // return;
    let info_div = d3.select("#pos_info_div").append("div")
    .style("text-align", "right")
    .style("padding-right", "40%");
    let categories = ["Active Site", "Ligand Binding", "Interface (Excluding Non-bonded)", "DNA Binding", "Disulfide bridge", "Salt bridge", "Co-Evo. Set"];
    let categories_list = [active_site_for_chain, all_ligand_binding_for_chain, interfacing_for_chain, dna_binding_for_chain, intrachain_disulphide_list, salt_bridges, coevolution_pos];
    for (let i = 0; i < categories.length; i++) {
        const category_name = categories[i];
        const category_list = categories_list[i];
        let is_in_category = "✗";
        let category_color = "red";
        if (category_list.indexOf(trespos) > -1) {
            is_in_category = "✔";
            category_color = "green";
        }
        if (category_name === "Salt bridge" && category_list.indexOf(resp_num+'') > -1) {
            is_in_category = "✔";
            category_color = "green";
        }
        let category_to_table = info_div.append("span")
        .style("font-family", "sans-serif")
        .style("font-size", "20px")
        ;
        category_to_table.append("span")
        .html(`${category_name}`);

        category_to_table.append("span")
        .style("color", category_color)
        .style("padding-left", "10%")
        .html(`${is_in_category}`)
        info_div.append("br");
    }
}

function generate_mutation_view(mutpos) {
    let d3_mutation_view = d3.select("#mutation_view");
    d3_mutation_view.html("");
    let content = d3_mutation_view.html();
    if (!content) {
        d3_mutation_view.style("text-align", "center");
        d3_mutation_view.style("max-height", width_80_perc+"px");
        d3_mutation_view.style("margin-left", width_10_perc+"px");
        d3_mutation_view.style("margin-right", width_10_perc+"px");
        d3_mutation_view.append("br");
        d3_mutation_view.append("br");
        d3_mutation_view.append("hr");
        d3_mutation_view.append("h1")
        .html("Mutation Info")
        ;
        d3_mutation_view.append("div")
        .attr("id", "mut_info_div");
    } else {
        console.log("exists");
        d3.select("#mut_info_div").html("");
    }
    let mutpos_num = mutpos.substr(1, mutpos.length-2);
    let wt_res = mutpos[0]
    let mut_res = mutpos[mutpos.length-1];
    let twt_res = oneToThree[wt_res];
    let tmut_res = oneToThree[mut_res];
    let pdb_mut = twt_res + mutpos_num + tmut_res;
    let unipresp_num = uniprot_residue_mapping[twt_res+mutpos_num];
    let up_mut = twt_res + unipresp_num + tmut_res;
    console.log("rectangles_array");
    console.log(rectangles_array[0]);
    console.log(rectangles_array[1]);
    console.log(rectangles_array[2]);
    let rectangles_array_f = rectangles_array.filter(function(a) {
        return a.c !== "None";
    });
    let mutations_rect_sorted = rectangles_array_f.sort(function(a, b) {
            return b.c - a.c;
    });
    console.log("mutations_rect_sorted");
    console.log(mutations_rect_sorted[0]);
    console.log(mutations_rect_sorted[1]);
    console.log(mutations_rect_sorted[2]);
    let mutations_rect_ids = mutations_rect_sorted.map(function(a) {return a.mutpos});
    let mut_rank = mutations_rect_ids.indexOf(mutpos);

    d3.select("#mut_info_div").append("h2")
    .html(`Mutation (PDB): ${pdb_mut}`)
    ;

    d3.select("#mut_info_div").append("h2")
    .html(`Mutation(UniProt): ${up_mut}`)
    ;

    var shap_data = prediction_per_mut[mutpos];
    if (shap_data[1] === "None") {
        return;
    }
    var shap_detailed = shap_data[2];

    d3.select("#mut_info_div").append("h2")
    .html(`Score: ${formatVal(shap_data[1], 3)}`)
    ;

    d3.select("#mut_info_div").append("h2")
    .html(`Score Rank: ${mut_rank+1} of ${mutations_rect_ids.length}`)
    .append("span")
    .html(' (Higher → More deleterious)').style("font-size", "16px");
    ;
    d3.select("#mut_info_div").append("br");
    d3.select("#mut_info_div").append("br");

    d3.select("#mut_info_div").append("h2")
    .html("SHAP Feature Importance:")
    ;

    d3.select("#mut_info_div").append("h3")
    .html(`Base Prediction value: ${formatVal(shap_detailed.base, 3)}`)
    ;

    d3.select("#mut_info_div").append("h3")
    .html(`Contributions of features: ${formatVal((shap_data[1]-shap_detailed.base), 3)}`)
    ;

    console.log("shap_data");
    console.log(shap_data);
    console.log("shap_detailed");
    console.log(shap_detailed);

    var shap_plot_svg_axes = d3.select("#mut_info_div").append("svg")
    .attr("id", "svg_shapplot_axes")
    .attr("width", total_width_general_s/2)
    .attr("height", "300px");

    var shap_plot_svg = d3.select("#mut_info_div").append("svg")
    .attr("id", "svg_shapplot")
    .attr("width", total_width_general_s+line_plot_sx*2)
    .attr("height", "300px");

    // var shap_detailed_values = shap_detailed.values.map(function(a) {
    //     // return -1.0 * a;
    //     return 1.0 * a;
    // });
    var shap_detailed_values = shap_detailed.values

    let shap_x_extent = d3.extent(shap_detailed_values);
    let shap_x_extent2 = [Math.abs(shap_x_extent[0]),Math.abs(shap_x_extent[1])];
    console.log(shap_x_extent);
    console.log(shap_x_extent2);
    console.log(d3.max(shap_x_extent2));
    // shap_x_extent2 = shap_x_extent2.sort(function(a,b) {
    //     return b-a;
    // });
    let half_shap_plot = total_width_general_s/2;
    var plot_remainder = total_width_general_s-(half_shap_plot+line_plot_sx*2);

    var shapScaleX = d3.scaleLinear()
    .domain([0.0, d3.max(shap_x_extent2)])
    .range([0, half_shap_plot-line_plot_sx]);

    // var shapScaleX;
    // if (shap_x_extent[0] < 0.0  && shap_x_extent[1] > 0.0) {
    //     console.log("1 scale");
    //     shapScaleX = d3.scaleLinear()
    //     .domain([shap_x_extent[0],0.0,shap_x_extent[1]])
    //     .range([line_plot_sx, half_shap_plot, total_width_general_s-line_plot_sx]);
    // } else if (shap_x_extent[0] < 0.0  && shap_x_extent[1] < 0.0) {
    //     console.log("2 scale");
    //     shapScaleX = d3.scaleLinear()
    //     .domain([shap_x_extent[1], shap_x_extent[0]])
    //     .range([line_plot_sx,total_width_general_s-line_plot_sx]);
    // } else if (shap_x_extent[0] > 0.0  && shap_x_extent[1] > 0.0) {
    //     console.log("3 scale");
    //     shapScaleX = d3.scaleLinear()
    //     .domain([shap_x_extent[0], shap_x_extent[1]])
    //     .range([line_plot_sx,total_width_general_s-line_plot_sx]);
    // }
    // console.log(shapScaleX(shap_x_extent[0]));
    // console.log(shapScaleX(0.0));
    // console.log(shapScaleX(shap_x_extent[1]));

    var shapScaleY = d3.scaleLinear()
    .domain([0,shap_detailed_values.length])
    .range([20,280]);

    shap_plot_svg.selectAll(".shap_bars")
    .data(shap_detailed_values)
    .enter()
    .append("rect")
    .attr("class", "shap_bars")
    // .attr("x",function(d,i) {return shapScaleX(0.0)})
    .attr("x",function(d,i) {return half_shap_plot})
    .attr("y", function(d,i) {return shapScaleY(i)})
    .attr("width", function(d) {
        if (d<0) {
            return shapScaleX(-1.0*d);
        }
        return shapScaleX(d);
    })
    .attr("height", function() {return shapScaleY(0)})
    .attr("fill", function(d) {
        if (d > 0) {
            return "red";
        } else {
            return "blue";
        }
    })
    .style("transform-box", "fill-box")
    .style("transform-origin", "center left")
    .style("transform", function(d) {
        if (d < 0) {
            return "rotate(180deg)"
        }
        return;
    })
    ;

    shap_plot_svg.selectAll(".shap_bars_text")
    .data(shap_detailed_values)
    .enter()
    .append("text")
    .attr("class", "shap_bars_text")
    .attr("x",function(d,i) {return half_shap_plot})
    .attr("y", function(d,i) {return shapScaleY(i)})
    .text(function(d) {
        return formatVal(d, 3);
    });

    shap_plot_svg_axes.selectAll(".shap_bars_text")
    .data(shap_detailed.names)
    .enter()
    .append("text")
    .attr("class", "shap_bars_text")
    .attr("x", total_width_general_s/6)
    .attr("y", function(d,i) {return shapScaleY(i+0.5)})
    .text(function(d) {
        return d;
    });
}
