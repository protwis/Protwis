// #################
// ###   TREE    ###
// #################

// Restructure tree data if only one Class is present
function update_tree_data(data,depth) {
    if (depth === 4) {
        // Iterate over each class and update name
        data.children.forEach(Class_child => {
            Class_child.name = Class_child.name.split(" (")[0];
            // iterate over ligand type
            Class_child.children.forEach(LigandType => {
                // iterate over receptor family
                LigandType.children.forEach(ReceptorFamily => {
                    // Trimming receptor family
                    ReceptorFamily.name = ReceptorFamily.name.replace(/( receptors|neuropeptide )/g, '');
                    ReceptorFamily.name = ReceptorFamily.name.split(" (")[0];
                });
            });

        });
    } else if (depth === 3) {
         // iterate over ligand type
         data.children.forEach(LigandType => {
            // iterate over receptor family
            LigandType.children.forEach(ReceptorFamily => {
                // Trimming receptor family
                ReceptorFamily.name = ReceptorFamily.name.replace(/( receptors|neuropeptide )/g, '');
                ReceptorFamily.name = ReceptorFamily.name.split(" (")[0];
            });
        });
    }
    return data
}

// Draw the Tree
function draw_tree(data, options,circle_size) {

    // Remove existing SVG if present
    d3.select('#' + options.anchor + "_svg").remove();

    var branches = {};
    var branch_offset = 0;
    var thickness = options.depth + 1;
    for (var key in options.branch_length) {
        if (key === options.depth) { continue };
        if (options.label_free.includes(parseInt(key,10))) {
            branch_offset = branch_offset + 10;
        } else {
            if (options.branch_trunc !== 0) {
                branch_offset = branch_offset + 2 * options.branch_trunc + 10;
            } else {
                branch_offset = branch_offset + string_pixlen(options.branch_length[key], key);
            }
        }
        branches[key] = branch_offset;
    }
    branches[options.depth] = branch_offset + options.leaf_offset;

    var diameter = 2 * branches[options.depth] + 140;

    var tree = d3.layout.tree()
        .size([360, diameter / 2])
        .separation(function (a, b) { return (a.parent === b.parent ? 1 : 2) / a.depth; });

    var diagonal = d3.svg.diagonal.radial()
        .projection(function (d) { return [d.y, d.x / 180 * Math.PI]; });

    var svg = d3.select('#' + options.anchor).append("svg")
        .attr("width", diameter)
        .attr("height", diameter)
        .attr("id", options.anchor + "_svg")
        .attr("xmlns", "http://www.w3.org/2000/svg");

    var svg_g = svg.append("g")
        .attr("transform", "translate(" + diameter / 2 + "," + diameter / 2 + ")");

    var nodes = tree.nodes(data);

    nodes.forEach(function (d) {
        if (d.depth === 0) {
            d.y = 0
        } else {
            d.y = branches[d.depth]
        }
    });

    var links = tree.links(nodes);

    var link = svg_g.append("g")
        .attr("class", "links")
        .selectAll("path")
        .data(links)
        .enter().append("path")
        .each(function (d) { d.target.linkNode = this; })
        .attr("d", diagonal) //function (d) { return step(d.source.x, d.source.y, d.target.x, d.target.y) })
        .style("stroke", function (d) { return d.target.color; })
        .style("stroke-width", function (d) { if (d.target.depth > 0) { return thickness - d.target.depth; } else { return 0; } })
        .style("fill-opacity", 0)
        .style("opacity", function (d) {
            if ((d.target.interactions > 0 && d.target.mutations_an > 0) || 1 === 1) { return 0.8 } //|| 1==1
            else if (d.target.interactions > 0) { return 0.5 }
            else if (d.target.mutations_an > 0) { return 0.5 }
            else { return 0.1 };
        });

    var node = svg_g.selectAll(".node")
        .data(nodes)
        .enter().append("g")
        .attr("class", "node")
        .attr("transform", function (d) { if (d.name === '') { return "rotate(" + (d.x) + ")translate(" + d.y + ")"; } else { return "rotate(" + (d.x - 90) + ")translate(" + d.y + ")"; } })

    node.filter(function (d) { return (d.depth === options.depth) })
        .attr("id", function (d) { if (d.name === '') { return "innerNode" } else { return 'X' + d.name.toUpperCase() } });

    node.append("text")
        .attr("dy", ".31em")
        .attr("name", function (d) { if (d.name === '') { return "branch" } else { return d.name } })
        .attr("text-anchor", function (d) {
            if (d.depth === options.depth) {
                return d.x < 180 ? "start" : "end";
            } else {
                return d.x < 180 ? "end" : "start";
            }
        })
        .attr("transform", function (d) {
            var labelOffset = parseFloat(circle_size) + 4;  // Adjust the offset as needed
            if (d.depth === options.depth) {
                return d.x < 180 ? `translate(${labelOffset})` : `rotate(180)translate(-${labelOffset})`;
            } else {
                return d.x < 180 ? "translate(-12)" : "rotate(180)translate(12)";
            }
        })
        .text(function (d) {
            if (d.depth === options.depth) {
                return d.name.toUpperCase();
            } else if (options.label_free.includes(d.depth)) {
                return "";
            } else if (d.depth > 0) {
                return d.name;
            } else {
                return "";
            }
        })
        .call(wrap, options.branch_trunc)
        .style("font-size", function (d) {
            // Use the custom font size from options
            if (options.depth === 4) {
            if (d.depth === 1) { return options.fontSize.class; }
            else if (d.depth === 2) { return options.fontSize.ligandtype; }
            else if (d.depth === 3) { return options.fontSize.receptorfamily; }
            else { return options.fontSize.receptor; }
        } else {
            if (d.depth === 1) { return options.fontSize.ligandtype; }
            else if (d.depth === 2) { return options.fontSize.receptorfamily; }
            else { return options.fontSize.receptor; }
        }
        })
        .style("font-family", "Palatino")
        .style("fill", function (d) {
            if (d.color) { return "#111"; }
            else { return "#222"; };
        }).call(getBB);

    node.filter(function (d) { return (d.depth !== options.depth) }).insert("rect", "text")
        .attr("x", function (d) { return d.x < 180 ? d.bbox.x - 12 : d.bbox.x - d.bbox.width - 12; })
        .attr("y", function (d) { return d.bbox.y; })
        .attr("width", function (d) { return d.bbox.width; })
        .attr("height", function (d) { return d.bbox.height; })
        .style("fill", "#FFF");

    function step(startAngle, startRadius, endAngle, endRadius) {
        var c0 = Math.cos(startAngle = (startAngle - 90) / 180 * Math.PI),
            s0 = Math.sin(startAngle),
            c1 = Math.cos(endAngle = (endAngle - 90) / 180 * Math.PI),
            s1 = Math.sin(endAngle);
        return "M" + startRadius * c0 + "," + startRadius * s0
            + (endAngle === startAngle ? "" : "A" + startRadius + "," + startRadius + " 0 0 " + (endAngle > startAngle ? 1 : 0) + " " + startRadius * c1 + "," + startRadius * s1)
            + "L" + endRadius * c1 + "," + endRadius * s1;
    }

    function string_pixlen(text, depth) {
        var canvas = document.createElement('canvas');
        var ctx = canvas.getContext("2d");

        // Check the depth condition
        if (options.depth === 4) {
            // Use the custom font size from options for all levels
            if (depth === 1) {
                ctx.font = options.fontSize.class + " Palatino";
            } else if (depth === 2) {
                ctx.font = options.fontSize.ligandtype + " Palatino";
            } else if (depth === 3) {
                ctx.font = options.fontSize.receptorfamily + " Palatino";
            } else {
                ctx.font = options.fontSize.receptor + " Palatino";
            }
        } else if (options.depth === 3) {
            // Omit class and start from ligandtype
            if (depth === 1) {
                ctx.font = options.fontSize.ligandtype + " Palatino";
            } else if (depth === 2) {
                ctx.font = options.fontSize.receptorfamily + " Palatino";
            } else {
                ctx.font = options.fontSize.receptor + " Palatino";
            }
        }

        // Measure and return the text width plus padding
        return parseInt(ctx.measureText(text).width,10);
    }

    function getBB(selection) {
        selection.each(function (d) { d.bbox = this.getBBox(); })
    }

    function wrap(text, width) {
        if (width === 0) {
            return;
        }
        text.each(function () {
            var text = d3.select(this),
                words = text.text().split(/\s+/).reverse(),
                word,
                line = [],
                lineNumber = 0,
                lineHeight = 1.1, // ems
                y = text.attr("y"),
                dy = parseFloat(text.attr("dy")),
                tspan = text.text(null).append("tspan").attr("x", 0).attr("y", y).attr("dy", dy + "em");

            word = words.pop(); // Initialize the first word
            while (word !== undefined) {
                line.push(word);
                tspan.text(line.join(" "));

                if (tspan.node().getComputedTextLength() > width) {
                    line.pop(); // Remove the last word as it exceeds the width
                    tspan.text(line.join(" ")); // Update the tspan with the valid words
                    line = [word]; // Start a new line with the word that didn't fit
                    tspan = text.append("tspan")
                                .attr("x", 0)
                                .attr("y", y)
                                .attr("dy", ++lineNumber * lineHeight + dy + "em")
                                .text(word); // Add the word to the new tspan
                }

                word = words.pop(); // Reassign the next word at the end of the loop
            }
        });
    }
    // === Centering Logic ===
    var scaleFactor = 0.8;  // Adjust this as needed

    // Calculate the extra padding needed based on circle size and spacer
    var extraPadding = (circle_size) + (circle_spacer*2) + parseInt(options.fontSize.receptor,10)*4;  // Adjust the multiplier based on how much padding is needed

    // Calculate new dimensions
    var newWidth = diameter + extraPadding;  // Add padding to both sides
    var newHeight = diameter + extraPadding;  // Add padding to both sides

    // Set new width and height for the SVG
    svg.attr('width', newWidth)
       .attr('height', newHeight);

    // Calculate the center point
    var cx = newWidth / 2;
    var cy = newHeight / 2;

    // Adjust the translation to keep the tree centered after scaling
    var translateX = cx;
    var translateY = cy;

    // Apply the transform to the 'g' element to scale and center it
    svg.select('g')
        .attr('transform', `translate(${translateX},${translateY}) scale(${scaleFactor},${scaleFactor})`);
}

// replaces labels derived from view
function formatTextWithHTML(text) {
    // Apply all the replacements step by step
    return text
        .replace(" receptor", '')
        .replace("-adrenoceptor", '')
        .replace(" receptor-", '-')
        .replace("<sub>", '</tspan><tspan baseline-shift="sub">')
        .replace("</sub>", '</tspan><tspan>')
        .replace("<i>", '</tspan><tspan font-style="italic">')
        .replace("</i>", '</tspan><tspan>')
        .replace("Long-wave-sensitive", 'LWS')
        .replace("Medium-wave-sensitive", 'MWS')
        .replace("Short-wave-sensitive", 'SWS')
        .replace("Olfactory", 'OLF')
        .replace("calcitonin-like receptor", 'CLR');
}

// Change the labels
function changeLeavesLabels(location, value, dict){
    // Initialize leaf node length
    maxLeafNodeLenght = 0;
    // Find longest label
    gNodes = d3.select('#'+location).selectAll('g');
    gNodes.each(function(d) {
      if (d3.select(this).attr("id") !== null) {
        name = d3.select(this).attr("id").substring(1);
        labelName = dict[name][0];
        labelName = formatTextWithHTML(labelName)
        node = d3.select('#'+location).select('#X'+name);
        if (node.size() !== 0){
          if (value === "IUPHAR"){
            node.selectAll("text")[0].forEach(
              function(node_label){
                node_label.innerHTML = labelName;
                labelSize = node_label.getBBox().width*1.05 + 0.5 * 10
                if (labelSize > maxLeafNodeLenght){
                  // change initialization label length, needed for outer circles
                  maxLeafNodeLenght = labelSize
                }
              });
          } else if (value === "UniProt"){
            node.selectAll("text")[0].forEach(
              function(node_label){
                node_label.innerHTML = name;
                labelSize = node_label.getBBox().width*1.05 + 0.5 * 10
                if (labelSize > maxLeafNodeLenght){
                  maxLeafNodeLenght = labelSize
                }
              });
          }
        }
      }
    });
  }

// Draw the circles (data) of the tree plot

function DrawCircles(location, data, starter, dict, clean = true, gradient = true, circle_styling_dict, circle_spacer, circle_size) {
    var svg = d3.select('#' + location);
    var node = svg.selectAll(".node");

    if (clean === true) {
        node.selectAll("circle.outerCircle, circle.innerCircle").remove();  // Remove previously drawn circles (both outer and inner)
    }

    var spacer = circle_spacer;

    // Initialize a dictionary to store min and max values for each unit
    var minMaxValues = {};

    // First pass: Determine min and max values for each unit (including 'Inner')
    for (var x in data) {
        var keys = Object.keys(data[x]);
        for (var unit of keys) {
            if (!minMaxValues[unit]) {
                minMaxValues[unit] = { min: Infinity, max: -Infinity };
            }
            var value = data[x][unit];
            if (value < minMaxValues[unit].min) {
                minMaxValues[unit].min = value;
            }
            if (value > minMaxValues[unit].max) {
                minMaxValues[unit].max = value;
            }
        }
    }

    // Second pass: Draw the circles for both "Inner" and "Outer" units
    for (var x in data) {
        var keys = Object.keys(data[x]);

        for (var unit of keys) {
            var leaf = svg.selectAll('g[id=X' + x + ']');  // Use the node positions from the tree

            if (unit === 'Inner') {
                // Draw 'Inner' circle at the node position
                if (dict[unit]) {
                    var value = data[x][unit];
                    var minValue = minMaxValues[unit].min;
                    var maxValue = minMaxValues[unit].max;

                    // Determine the styling for the 'Inner' unit
                    var styling = circle_styling_dict[unit] || "Two";

                    // Create color scale based on min and max values
                    var colorScale;
                    if (styling === "Three") {
                        // Three-color gradient with white in the middle
                        colorScale = d3.scale.linear()
                            .domain([minValue, (minValue + maxValue) / 2, maxValue])
                            .range([dict[unit][0], "#FFFFFF", dict[unit][1]]);
                    } else {
                        // Two-color gradient
                        colorScale = d3.scale.linear()
                            .domain([minValue, maxValue])
                            .range(dict[unit]);
                    }

                    // Calculate color using the color scale
                    var color = gradient ? colorScale(value) : dict[unit][0];  // Use the first color if no gradient

                    // Append the 'Inner' circle
                    leaf.append("circle")
                        .attr("r", circle_size)  // Draw 'Inner' circle
                        .attr("class", "innerCircle")  // Add class to distinguish inner circles
                        .style("stroke", "black")  // Use the first color for the stroke
                        .style("stroke-width", 0.8)
                        .style("fill", color)
                        .attr("transform", "translate(0,0)");  // Draw at the center of the node
                }
            } else {
                // Draw 'Outer' circles based on the unit
                if (dict[unit]) {
                    var value = data[x][unit];
                    var minValue = minMaxValues[unit].min;
                    var maxValue = minMaxValues[unit].max;

                    // Determine the styling for the outer units
                    var styling = circle_styling_dict[unit] || "Two";

                    // Create color scale based on min and max values
                    var colorScale;
                    if (styling === "Three") {
                        // Three-color gradient with white in the middle
                        colorScale = d3.scale.linear()
                            .domain([minValue, (minValue + maxValue) / 2, maxValue])
                            .range([dict[unit][0], "#FFFFFF", dict[unit][1]]);
                    } else {
                        // Two-color gradient
                        colorScale = d3.scale.linear()
                            .domain([minValue, maxValue])
                            .range(dict[unit]);
                    }

                    // Calculate color using the color scale
                    var color = gradient ? colorScale(value) : dict[unit][0];  // Use the first color if no gradient

                    var multiply = 1 + Object.keys(dict).indexOf(unit);  // For spacing between outer circles

                    // Append the outer circles, position them with `circle_spacer`
                    leaf.append("circle")
                        .attr("r", circle_size)
                        .attr("class", "outerCircle")  // Add class to distinguish outer circles
                        .style("stroke", "black")  // Use the first color for the stroke
                        .style("stroke-width", 0.8)
                        .style("fill", color)
                        .attr("transform", "translate(" + (Math.ceil(starter) + multiply * spacer) + ",0)");
                }
            }
        }
    }
}

// Create the bar legends
function createLegendBars(location, data, conversion, circle_styling_dict, datatype_dict) {
    var svg = d3.select('#' + location + ' svg');

    // Clear existing content
    svg.selectAll("*").remove();

    var margin = { top: 20, right: 20, bottom: 30, left: 60 };
    var width = +svg.attr("width") - margin.left - margin.right;
    var height = +svg.attr("height") - margin.top - margin.bottom;

    // Flatten data to get categories and their max values
    var categoryMax = {};

    Object.values(data).forEach(item => {
        Object.entries(item).forEach(([category, value]) => {
            if (category in conversion) {
                if (!categoryMax[category]) {
                    categoryMax[category] = { min: value, max: value };
                } else {
                    if (value > categoryMax[category].max) {
                        categoryMax[category].max = value;
                    }
                    if (value < categoryMax[category].min) {
                        categoryMax[category].min = value;
                    }
                }
            }
        });
    });

    var existingCategories = Object.keys(categoryMax).filter(cat => categoryMax[cat].max > 0);

    // Create a color scale function for each category
    var colorScales = {};
    existingCategories.forEach(category => {
        // Check the datatype in the datatype_dict
        if (datatype_dict[category] === "Discrete") {
            return; // Skip this category if its datatype is Discrete
        }

        var colors = conversion[category];
        var gradientId = `gradient-${category}`;

        var gradient = svg.append("defs")
            .append("linearGradient")
            .attr("id", gradientId)
            .attr("x1", "0%")
            .attr("x2", "100%")
            .attr("y1", "0%")
            .attr("y2", "0%");

        var styling = circle_styling_dict[category] || "Two";
        if (styling === "Three") {
            gradient.append("stop")
                .attr("offset", "0%")
                .attr("stop-color", colors[0]);
            gradient.append("stop")
                .attr("offset", "50%")
                .attr("stop-color", "#FFFFFF");
            gradient.append("stop")
                .attr("offset", "100%")
                .attr("stop-color", colors[1]);
        } else {
            gradient.append("stop")
                .attr("offset", "0%")
                .attr("stop-color", colors[0]);
            gradient.append("stop")
                .attr("offset", "100%")
                .attr("stop-color", colors[1]);
        }

        colorScales[category] = gradientId;
    });

    var barWidth = 120;
    var barHeight = 15;
    var spacing = 50; // Horizontal spacing between bars

    // Adjust SVG width if needed
    var totalWidth = 1200;
    svg.attr("width", totalWidth);

    var skippedDiscreteCount = 0;  // Track how many discrete bars have been skipped

    existingCategories.forEach((category, index) => {
        // Check if datatype is "Continuous" before drawing the legend
        if (datatype_dict[category] !== "Continuous") {
            skippedDiscreteCount++; // Increment the count of skipped discrete bars
            return; // Skip this category if it's Discrete
        }

        var gradientId = colorScales[category];
        var minValue = categoryMax[category].min;
        var maxValue = categoryMax[category].max;
        var midValue = (minValue + maxValue) / 2;
        minValue = parseFloat(minValue.toFixed(2));
        maxValue = parseFloat(maxValue.toFixed(2));
        midValue = parseFloat(midValue.toFixed(2));

        // Adjust the x position based on how many "Discrete" bars were skipped
        var xPosition = margin.left + (index - skippedDiscreteCount) * (barWidth + spacing);

        // Add a rectangle to represent the gradient
        svg.append("rect")
            .attr("x", xPosition)
            .attr("y", margin.top)
            .attr("width", barWidth)
            .attr("height", barHeight)
            .style("fill", `url(#${gradientId})`)
            .style("stroke", "black")
            .style("stroke-width", "1px");

        // Add a label for the gradient bar
        svg.append("text")
            .attr("x", xPosition + barWidth / 2)
            .attr("y", margin.top - 10)
            .attr("text-anchor", "middle")
            .text(category);

        // Add min and max value labels
        svg.append("text")
            .attr("x", xPosition)
            .attr("y", margin.top + barHeight + 15)
            .attr("text-anchor", "start")
            .text(`${minValue}`);

        svg.append("text")
            .attr("x", xPosition + barWidth)
            .attr("y", margin.top + barHeight + 15)
            .attr("text-anchor", "end")
            .text(`${maxValue}`);

        // Add mid value label
        svg.append("text")
            .attr("x", xPosition + barWidth / 2)
            .attr("y", margin.top + barHeight + 15)
            .attr("text-anchor", "middle")
            .text(`${midValue}`);
    });
}



// #################
// ###  CLUSTER  ###
// #################

// Function to naturally sort an array
function naturalSort(a, b) {
    return a.localeCompare(b, undefined, { numeric: true, sensitivity: 'base' });
}


// Function to create traces for the plot
function createTraces(colorOption, showLabels, colorMapping, textColorEnabled) {
    let traces = [];
    const stroke_width = cluster_DataStyling.strokeWidth;  // Get stroke width from data styling
    const marker_size = cluster_DataStyling.markerSize;    // Get marker size from data styling
    const border_on = cluster_DataStyling.borderOn;        // Get border toggle state from data styling

    const clusterSet = Array.from(new Set(currentClusterData.map(d => d.cluster))).sort((a, b) => a - b);

    // Handle clustering color option
    if (colorOption === 'cluster') {
        clusterSet.forEach(cluster => {
            const clusterData = currentClusterData.filter(d => d.cluster === cluster);

            // Create marker traces, only display markers when labels are not shown
            const markerTrace = {
                x: clusterData.map(d => d.x),
                y: clusterData.map(d => d.y),
                mode: showLabels ? 'none' : 'markers',  // Hide markers if labels are shown
                type: 'scatter',
                hoverinfo: 'text',
                text: clusterData.map(d => `${d.label.toUpperCase()}`),  // Combine label and fill for hover text
                marker: {
                    size: marker_size,
                    symbol: 'circle',
                    line: {
                        width: border_on ? stroke_width : 0,
                        color: 'black'
                    },
                    color: clusterData.map(d => colorPalette[d.cluster % colorPalette.length]),  // Keep legend colors based on clusters
                },
                // Customize legend text color only when showLabels and textColorEnabled are both true
                name: (showLabels && textColorEnabled)
                ? `<span style="color:${colorPalette[cluster % colorPalette.length]}">Cluster ${cluster + 1}</span>`
                : `Cluster ${cluster + 1}`  // Regular label if conditions are false
            };

            traces.push(markerTrace);
        });
    }
    // Handle gradient color option
    else if (colorOption === 'gradient') {
        const fillValues = currentClusterData.map(d => d.fill);
        const minFill = Math.min(...fillValues);
        const maxFill = Math.max(...fillValues);
        const gradientTrace = {
            x: currentClusterData.map(d => d.x),
            y: currentClusterData.map(d => d.y),
            mode: showLabels ? 'none' : 'markers',  // Hide markers if labels are shown
            type: 'scatter',
            hoverinfo: 'text',
            text: currentClusterData.map(d => `${d.label.toUpperCase()}: ${d.fill}`),  // Combine label and fill for hover text
            marker: {
                size: marker_size,
                symbol: 'circle',
                line: {
                    width: border_on ? stroke_width : 0,
                    color: 'black'
                },
                color: currentClusterData.map(d => d.fill),
                cmin: minFill,  // Set color axis minimum value
                cmax: maxFill,  // Set color axis maximum value
                colorscale: 'RdBu',  // Use RdBu color scale
                colorbar: {
                    thickness: 20,
                    len: 0.5,
                }
            },
            showlegend: false
        };

        traces.push(gradientTrace);
    }
    // Handle Class, Ligand type, or Receptor family color option
    else if (['Class', 'Ligand type', 'Receptor family'].includes(colorOption)) {

        const uniqueEntries = new Set();

        currentClusterData.forEach(point => {
            const entry = point[colorOption] === 'Other GPCRs' ? 'Classless' : point[colorOption];
            if (!uniqueEntries.has(entry)) {
                uniqueEntries.add(entry);

                // Create marker traces, only display markers when labels are not shown
                const markerTrace = {
                    x: currentClusterData.filter(d => d[colorOption] === point[colorOption]).map(d => d.x),
                    y: currentClusterData.filter(d => d[colorOption] === point[colorOption]).map(d => d.y),
                    mode: showLabels ? 'none' : 'markers',
                    type: 'scatter',
                    hoverinfo: 'text',
                    text: currentClusterData.filter(d => d[colorOption] === point[colorOption]).map(d => d.label.toUpperCase()),
                    marker: {
                        size: marker_size,
                        symbol: 'circle',
                        line: {
                            width: border_on ? stroke_width : 0,
                            color: 'black'
                        },
                        color: colorMapping[entry],
                    },
                    // Customize legend text color only when showLabels and textColorEnabled are both true
                    name: (showLabels && textColorEnabled)
                        ? `<span style="color:${colorMapping[entry]}">${entry}</span>`
                        : entry  // Regular label if conditions are false
                };

                traces.push(markerTrace);
            }
        });

        traces.sort((a, b) => naturalSort(a.name, b.name));
    }

    return traces;
}


// Function to create annotations for labels with optional text coloring based on `textColorEnabled`
function createAnnotations(filteredData, colorOption, textColorEnabled, colorMapping) {
    const annotations = [];

    // Get global min and max values for `fill` from all data (currentClusterData), not just the filtered data
    const fillValues = currentClusterData.map(d => d.fill);  // Use all data, not just filtered
    const minFill = Math.min(...fillValues);
    const maxFill = Math.max(...fillValues);

    // Use d3.interpolateRdBu for the exact RdBu color scale in D3 v4
    const rdBuColorScale = d3v4.scaleSequential(d3v4.interpolateRdBu)
        .domain([maxFill, minFill]);  // Inverse the domain for red to blue coloring


    filteredData.forEach((d) => {
        let textColor;

        if (textColorEnabled) {
            if (colorOption === 'cluster') {
                textColor = colorPalette[d.cluster % colorPalette.length];
            } else if (colorOption === 'gradient') {
                textColor = rdBuColorScale(d.fill);  // Gradient-based coloring
            } else if (['Class', 'Ligand type', 'Receptor family'].includes(colorOption)) {
                const entry = d[colorOption] === 'Other GPCRs' ? 'Classless' : d[colorOption];
                textColor = colorMapping[entry];
            } else {
                textColor = 'black';
            }
        } else {
            textColor = 'black';
        }

        annotations.push({
            x: d.x,
            y: d.y,
            xref: 'x',
            yref: 'y',
            text: `${d.label.toUpperCase()}`,  // Combine label and fill for hover text
            showarrow: false,
            font: {
                family: 'Arial',
                size: cluster_DataStyling.labelFontSize,
                color: textColor
            },
            align: 'center',
            bgcolor: 'rgba(0,0,0,0)',
            borderwidth: 0
        });
    });

    return annotations;
}

// Function to update the plot with markers or labels
function updatePlotWithAnnotations() {
    const colorOption = getActiveColorOption();  // Get the active color option
    const showLabels = document.getElementById('show').classList.contains('active');
    const textColorEnabled = cluster_DataStyling.textColorEnabled;

    const plotElement = document.getElementById('plotContainer_cluster');
    const currentLayout = plotElement ? Plotly.d3.select('#plotContainer_cluster').node().layout : {};

    const xRange = (currentLayout && currentLayout.xaxis && currentLayout.xaxis.range) ? currentLayout.xaxis.range : [Math.min(...currentClusterData.map(d => d.x)), Math.max(...currentClusterData.map(d => d.x))];
    const yRange = (currentLayout && currentLayout.yaxis && currentLayout.yaxis.range) ? currentLayout.yaxis.range : [Math.min(...currentClusterData.map(d => d.y)), Math.max(...currentClusterData.map(d => d.y))];

    let colorMapping = {};
    if (['Class', 'Ligand type', 'Receptor family'].includes(colorOption)) {
        let uniqueValues = Array.from(new Set(currentClusterData.map(d => d[colorOption])));
        uniqueValues = uniqueValues.map(value => value === 'Other GPCRs' ? 'Classless' : value);
        uniqueValues.sort(naturalSort);
        uniqueValues.forEach((value, index) => {
            colorMapping[value] = colorPalette[index % colorPalette.length];
        });
    }

    // Generate the traces and pass the colorMapping
    const traces = createTraces(colorOption, showLabels, colorMapping, textColorEnabled);

    // Add the color bar for gradient only if we're displaying annotations (text labels)
    if (colorOption === 'gradient' && showLabels) {
        const fillValues = currentClusterData.map(d => d.fill);
        const minFill = Math.min(...fillValues);
        const maxFill = Math.max(...fillValues);

        const colorbarTrace = {
            z: [[minFill, maxFill], [minFill, maxFill]],  // Use actual min/max values for z
            x: [0, 1],
            y: [0, 1],
            type: 'heatmap',
            colorscale: 'RdBu',
            showscale: true,  // Only show the color bar when annotations are visible
            colorbar: {
                thickness: 20,
                len: 0.5,
                // Removed title from the colorbar
            },
            opacity: 0  // Make the heatmap itself transparent
        };

        traces.push(colorbarTrace);
    }

    // Filter the data for labels within the current zoom range
    const filteredData = currentClusterData.filter(d => {
        const inXRange = (d.x >= xRange[0] && d.x <= xRange[1]);
        const inYRange = (d.y >= yRange[0] && d.y <= yRange[1]);
        return inXRange && inYRange;
    });

    // Generate annotations based on filtered data, text coloring state, and shared colorMapping
    const annotations = showLabels ? createAnnotations(filteredData, colorOption, textColorEnabled, colorMapping) : [];

    // Define new layout with annotations
    const layout = {
        xaxis: {
            visible: false,
            showgrid: false,
            range: xRange,
            scaleanchor: 'y'
        },
        yaxis: {
            visible: false,
            showgrid: false,
            range: yRange
        },
        hovermode: 'closest',
        showlegend: true,
        annotations: annotations,  // Add annotations to the plot
        legend: {
            x: 1,
            xanchor: 'left',
            y: 0.5,
            orientation: 'v'
        },
        plot_bgcolor: '#FFFFFF',
        autosize: false,
        width: 1024,
        height: 700,
        margin: {
            l: 50,
            r: 374,
            t: 50,
            b: 50
        }
    };

    Plotly.react('plotContainer_cluster', traces, layout);
}


// #################
// ###   LIST    ###
// #################

// ## Initialization of global data for list ##
function initializeDataStyling(list_data_wow, data_types_list) {

    // Function to check if any of the specified keys are present in a given object
    function checkKeysPresent(obj, keys) {
        return keys.some(k => k in obj);
    }

    // Define the keys to check for each column
    var col1Keys = ["Value1", "Value2"];
    var col2Keys = ["Value3", "Value4"];
    var col3Keys = ["Value5", "Value6"];
    var col4Keys = ["Value7", "Value8"];


    var Data_styling = {
        Col1: {Data: "No",Datatype: 'None',Data_min: Infinity,Data_max: -Infinity, Data_coloring: 'All', Data_color1: '#000000', Data_color2: '#000000', data_color_complexity: 'One'},
        Col2: {Data: "No",Datatype: 'None',Data_min: Infinity,Data_max: -Infinity, Data_coloring: 'All', Data_color1: '#000000', Data_color2: '#000000', data_color_complexity: 'One'},
        Col3: {Data: "No",Datatype: 'None',Data_min: Infinity,Data_max: -Infinity, Data_coloring: 'All', Data_color1: '#000000', Data_color2: '#000000', data_color_complexity: 'One'},
        Col4: {Data: "No",Datatype: 'None',Data_min: Infinity,Data_max: -Infinity, Data_coloring: 'All', Data_color1: '#000000', Data_color2: '#000000', data_color_complexity: 'One'}
        };

    // Iterate through initialData and update Data_styling accordingly
    Object.keys(list_data_wow).forEach(function(key) {
        var currentObj = list_data_wow[key];

        // Check for Col1 if not already set to "Yes"
        if (Data_styling.Col1.Data !== "Yes" && checkKeysPresent(currentObj, col1Keys)) {
            Data_styling.Col1.Data = "Yes";
        }

        // Check for Col2 if not already set to "Yes"
        if (Data_styling.Col2.Data !== "Yes" && checkKeysPresent(currentObj, col2Keys)) {
            Data_styling.Col2.Data = "Yes";
        }

        // Check for Col3 if not already set to "Yes"
        if (Data_styling.Col3.Data !== "Yes" && checkKeysPresent(currentObj, col3Keys)) {
            Data_styling.Col3.Data = "Yes";
        }

        // Check for Col4 if not already set to "Yes"
        if (Data_styling.Col4.Data !== "Yes" && checkKeysPresent(currentObj, col4Keys)) {
            Data_styling.Col4.Data = "Yes";
        }

        // If both Col1 and Col2 are "Yes", no need to continue iterating
        if (Data_styling.Col1.Data === "Yes" && Data_styling.Col2.Data === "Yes" && Data_styling.Col3.Data === "Yes" && Data_styling.Col4.Data === "Yes") {
            return; // Exit forEach loop early
        }
    });


    // Check for Col1 has data
    // Update Datatype for each column based on data_types
    ['Col1', 'Col2', 'Col3', 'Col4'].forEach(col => {
        if (Data_styling[col].Data === "Yes") {
            if (data_types_list[col] === 'Discrete') {
                Data_styling[col].Datatype = 'Discrete';
            } else if (data_types_list[col] === 'Continuous') {
                Data_styling[col].Datatype = 'Continuous';
            }
        }
    });


    // Run through data and set min and max values

    function updateDataMinMax(column, key, currentObj) {
        if (Data_styling[column].Data === "Yes" && Data_styling[column].Datatype === 'Continuous' && currentObj.hasOwnProperty(key)) {
            if (currentObj[key] > Data_styling[column].Data_max) {
                Data_styling[column].Data_max = currentObj[key];
            }
            if (currentObj[key] < Data_styling[column].Data_min) {
                Data_styling[column].Data_min = currentObj[key];
            }
        }
    }

    Object.keys(list_data_wow).forEach(function(key) {
        var currentObj = list_data_wow[key];

        updateDataMinMax('Col1', 'Value2', currentObj);
        updateDataMinMax('Col2', 'Value4', currentObj);
        updateDataMinMax('Col3', 'Value6', currentObj);
        updateDataMinMax('Col4', 'Value8', currentObj);
    });

    return Data_styling;
}

// Function to modify the data
function Initialize_Data(data) {
    // 1. Change "Other GPCRs" to "Classless" (Layer1)
    if (data["Other GPCRs"]) {
      data["Classless"] = data["Other GPCRs"];
      delete data["Other GPCRs"];
    }

    // 2. Change "Other GPCR orphans" to "Classless orphans" (Layer3)
    Object.keys(data).forEach(layer1Key => {
      const layer2Data = data[layer1Key];
      Object.keys(layer2Data).forEach(layer2Key => {
        const layer3Data = layer2Data[layer2Key];
        if (layer3Data["Other GPCR orphans"]) {
          layer3Data["Classless orphans"] = layer3Data["Other GPCR orphans"];
          delete layer3Data["Other GPCR orphans"];
        }
      });
    });

    // 3. Remove "Olfactory receptors" from Layer2
    Object.keys(data).forEach(layer1Key => {
      const layer2Data = data[layer1Key];
      if (layer2Data["Olfactory receptors"]) {
        delete layer2Data["Olfactory receptors"];
      }
    });

    return data;  // Return the modified data
  }

//   Create Classification dict
function Create_classification_dict(data) {
    // Initialize the arrays for present items
    let present_class = [];
    let present_LigandType = [];
    let present_ReceptorFamilies = [];
    let present_receptors = [];

    // Iterate over the nested dictionary to split into the respective arrays
    for (const classKey in data) {
        present_class.push(classKey);  // Add the class key to the present_class array

        for (const ligandKey in data[classKey]) {
            present_LigandType.push(ligandKey);  // Add the ligand type key to the present_LigandType array

            for (const receptorFamilyKey in data[classKey][ligandKey]) {
                present_ReceptorFamilies.push(receptorFamilyKey);  // Add receptor family key to present_ReceptorFamilies

                // Add the receptors to the present_receptors array
                present_receptors.push(...data[classKey][ligandKey][receptorFamilyKey]);
            }
        }
    }

}

// Create array of sorted entries and classification array for printing
function Data_resorter(data) {

    // Check the visibility of each layer
    const Layer1_isChecked = d3.select(`#toggle-layer-1`).property('checked');
    const Layer2_isChecked = d3.select(`#toggle-layer-2`).property('checked');
    const Layer3_isChecked = d3.select(`#toggle-layer-3`).property('checked');

    // Data sorter depending on layers being shown:

    // Initialize the final array for results and the category array for categories
    let final_array = [];
    let category_array = [];

    // Helper function to sort arrays naturally
    function naturalSort(arr) {
        return arr.sort((a, b) => a.localeCompare(b, undefined, { numeric: true }));
    }

    // Build the sorted Class list (only include it if Layer1 is checked)
    let sortedClasses = Layer1_isChecked ? naturalSort(Object.keys(data)) : [];

    // Iterate over sorted Classes and handle layers accordingly
    for (const classKey of sortedClasses) {
        final_array.push(classKey); // Add the Class if Layer1 is checked
        category_array.push("Class"); // Add "Class" category

        let ligandTypeDict = {};

        // Collect all LigandTypes and their families/receptors for the current Class
        for (const ligandKey in data[classKey]) {
            if (!ligandTypeDict[ligandKey]) {
                ligandTypeDict[ligandKey] = [];
            }

            for (const receptorFamilyKey in data[classKey][ligandKey]) {
                ligandTypeDict[ligandKey].push({
                    family: receptorFamilyKey,
                    receptors: data[classKey][ligandKey][receptorFamilyKey]
                });
            }
        }

        // Sort LigandTypes if Layer2 is checked, otherwise skip LigandTypes
        let sortedLigandTypes = Layer2_isChecked ? naturalSort(Object.keys(ligandTypeDict)) : [];

        // For each LigandType, process its ReceptorFamilies and Receptors
        sortedLigandTypes.forEach(ligandKey => {
            final_array.push(ligandKey); // Add LigandType if Layer2 is checked
            category_array.push("LigandType"); // Add "LigandType" category

            let receptorFamilies = ligandTypeDict[ligandKey].sort((a, b) => a.family.localeCompare(b.family)); // Sort ReceptorFamilies

            receptorFamilies.forEach(familyObj => {
                if (Layer3_isChecked) {
                    final_array.push(familyObj.family); // Add ReceptorFamily if Layer3 is checked
                    category_array.push("ReceptorFamily"); // Add "ReceptorFamily" category
                }

                // Always sort and add receptors
                familyObj.receptors.forEach(receptor => {
                    final_array.push(receptor);
                    category_array.push("Receptor"); // Add "Receptor" category for each receptor
                });
            });
        });

        // If Layer2 is not checked but Layer3 is checked (Case [True, False, True])
        if (!Layer2_isChecked && Layer3_isChecked) {
            let receptorFamilies = [];
            for (const ligandKey in ligandTypeDict) {
                receptorFamilies.push(...ligandTypeDict[ligandKey]);
            }

            receptorFamilies.sort((a, b) => a.family.localeCompare(b.family));

            receptorFamilies.forEach(familyObj => {
                final_array.push(familyObj.family); // Add ReceptorFamily if Layer3 is checked
                category_array.push("ReceptorFamily"); // Add "ReceptorFamily" category
                familyObj.receptors.forEach(receptor => {
                    final_array.push(receptor);
                    category_array.push("Receptor"); // Add "Receptor" category
                });
            });
        }

        if (!Layer2_isChecked && !Layer3_isChecked) {
            for (const ligandKey in ligandTypeDict) {
                ligandTypeDict[ligandKey].forEach(familyObj => {
                    familyObj.receptors.forEach(receptor => {
                        final_array.push(receptor); // Add sorted Receptors
                        category_array.push("Receptor"); // Add "Receptor" category
                    });
                });
            }
        }
    }

    // Case [False, True, True]: Include LigandType, ReceptorFamily, and Receptors, but no Class
    if (!Layer1_isChecked && Layer2_isChecked && Layer3_isChecked) {
        let allLigandTypes = [];

        // Collect LigandTypes and ReceptorFamilies across all classes
        for (const classKey in data) {
            for (const ligandKey in data[classKey]) {
                if (!allLigandTypes[ligandKey]) {
                    allLigandTypes[ligandKey] = [];
                }
                for (const receptorFamilyKey in data[classKey][ligandKey]) {
                    allLigandTypes[ligandKey].push({
                        family: receptorFamilyKey,
                        receptors: data[classKey][ligandKey][receptorFamilyKey]
                    });
                }
            }
        }

        let sortedLigandTypes = naturalSort(Object.keys(allLigandTypes));

        sortedLigandTypes.forEach(ligandKey => {
            final_array.push(ligandKey); // Add LigandType
            category_array.push("LigandType"); // Add "LigandType" category

            let receptorFamilies = allLigandTypes[ligandKey].sort((a, b) => a.family.localeCompare(b.family));

            receptorFamilies.forEach(familyObj => {
                final_array.push(familyObj.family); // Add ReceptorFamily
                category_array.push("ReceptorFamily"); // Add "ReceptorFamily" category

                familyObj.receptors.forEach(receptor => {
                    final_array.push(receptor);
                    category_array.push("Receptor"); // Add "Receptor" category
                });
            });
        });
    }

    // Case [False, True, False]: Include LigandType and Receptors, but no Class or ReceptorFamily
    if (!Layer1_isChecked && Layer2_isChecked && !Layer3_isChecked) {
        let allLigandTypes = [];

        // Collect LigandTypes and Receptors across all classes
        for (const classKey in data) {
            for (const ligandKey in data[classKey]) {
                if (!allLigandTypes[ligandKey]) {
                    allLigandTypes[ligandKey] = [];
                }
                for (const receptorFamilyKey in data[classKey][ligandKey]) {
                    allLigandTypes[ligandKey].push(...data[classKey][ligandKey][receptorFamilyKey]);
                }
            }
        }

        let sortedLigandTypes = naturalSort(Object.keys(allLigandTypes));

        sortedLigandTypes.forEach(ligandKey => {
            final_array.push(ligandKey); // Add LigandType
            category_array.push("LigandType"); // Add "LigandType" category

            naturalSort(allLigandTypes[ligandKey]).forEach(receptor => {
                final_array.push(receptor);
                category_array.push("Receptor"); // Add "Receptor" category
            });
        });
    }

    // If no Class or LigandType is checked but Receptors need to be included
    if (!Layer1_isChecked && !Layer2_isChecked && Layer3_isChecked) {
        let allReceptorFamilies = [];

        for (const classKey in data) {
            for (const ligandKey in data[classKey]) {
                for (const receptorFamilyKey in data[classKey][ligandKey]) {
                    // Collect receptor families and their receptors directly from data
                    allReceptorFamilies.push({
                        family: receptorFamilyKey,
                        receptors: data[classKey][ligandKey][receptorFamilyKey]
                    });
                }
            }
        }

        // Sort the receptor families
        allReceptorFamilies.sort((a, b) => a.family.localeCompare(b.family));

        allReceptorFamilies.forEach(familyObj => {
            final_array.push(familyObj.family); // Add ReceptorFamily
            category_array.push("ReceptorFamily"); // Add "ReceptorFamily" category
            familyObj.receptors.forEach(receptor => {
                final_array.push(receptor);
                category_array.push("Receptor"); // Add "Receptor" category
            });
        });
    }


    // If all layers are unchecked, add all sorted receptors
    if (!Layer1_isChecked && !Layer2_isChecked && !Layer3_isChecked) {
        for (const classKey in data) {
            for (const ligandKey in data[classKey]) {
                for (const receptorFamilyKey in data[classKey][ligandKey]) {
                    // Access receptors directly and sort them
                    naturalSort(data[classKey][ligandKey][receptorFamilyKey]).forEach(receptor => {
                        final_array.push(receptor); // Add sorted Receptors
                        category_array.push("Receptor"); // Add "Receptor" category
                    });
                }
            }
        }
    }


    return { final_array, category_array };
}

// Calculate the dimensions of the plot
function Calculate_dimension(data, Category_data, Col_break_number, columns, label_conversion_dict, label_names, styling_option) {

    // Define the column tracking variables
    let temp_col_state = 1; // Current column
    let label_dim_counter = 0; // Counter for how many labels processed in the current column

    // Initialize label max width tracking for up to 4 columns
    let label_max_dict = {
        col1_label_max: { 'Class': -Infinity, 'LigandType': -Infinity, 'ReceptorFamily': -Infinity, 'Receptor': -Infinity },
        col2_label_max: { 'Class': -Infinity, 'LigandType': -Infinity, 'ReceptorFamily': -Infinity, 'Receptor': -Infinity },
        col3_label_max: { 'Class': -Infinity, 'LigandType': -Infinity, 'ReceptorFamily': -Infinity, 'Receptor': -Infinity },
        col4_label_max: { 'Class': -Infinity, 'LigandType': -Infinity, 'ReceptorFamily': -Infinity, 'Receptor': -Infinity }
    };

    let Col_spacing_dict = { 'Col1': -Infinity, 'Col2': -Infinity, 'Col3': -Infinity, 'Col4': -Infinity };

    // List of known HTML entities to replace for correct label width calculation (for IUPHAR)
    const htmlEntities = ['&alpha;', '&beta;', '&gamma;', '&delta;', '&epsilon;', '&zeta;', '&eta;', '&theta;', '&iota;', '&kappa;', '&lambda;', '&mu;', '&nu;', '&xi;', '&omicron;', '&pi;', '&rho;', '&sigma;', '&tau;', '&upsilon;', '&phi;', '&chi;', '&psi;', '&omega;'];

    // Helper function to replace HTML entities in label text (only for IUPHAR)
    function replaceHtmlEntities(str) {
        const entityRegex = new RegExp(htmlEntities.join('|'), 'gi');
        return str.replace(entityRegex, 'a'); // Replace with 'a' or some neutral character of approximate size
    }

    // Process each label in the data and check its corresponding category from Category_data
    data.forEach((label, index) => {
        const category = Category_data[index]; // Get the category of the current label ("Class", "LigandType", "ReceptorFamily", or "Receptor")

        // Increment the label dimension counter
        label_dim_counter++;

        // Move to the next column if the Col_break_number is exceeded
        if (label_dim_counter > Col_break_number && temp_col_state < columns) {
            label_dim_counter = 1; // Reset the counter
            temp_col_state++; // Move to the next column
        }

        // Trimming and processing labels based on category
        if (category === 'Class') {
            // Trimming Class (e.g., "Class (some text)" becomes "Class")
            label = label.split(" (")[0];

        } else if (category === 'LigandType') {
            // No additional processing for LigandType
            label = label;

        } else if (category === 'ReceptorFamily') {
            // Trimming ReceptorFamily (removing "receptors" or "neuropeptide" and trimming after "(")
            label = label.replace(/( receptors|neuropeptide )/g, '').split(" (")[0];

        } else if (category === 'Receptor') {
            // For Receptors, apply the label conversion based on 'label_names'
            if (label_names === 'Gene') {
                // Convert based on UniProt data
                label = label_conversion_dict[label];
                label = label ? label.replace(/_human/g, '').toUpperCase() : label; // Clean up UniProt receptor names
            } else if (label_names === 'Protein') {
                // Apply IUPHAR-specific replacements and clean up
                label = replaceHtmlEntities(label)
                    .replace(/<\/?i>|(-adrenoceptor| receptor)|<\/?sub>/g, '');

                // Subscript handling for accurate label length measurement
                label = label.replace(/<sub>.*?<\/sub>/g, ''); // Simplified subscript removal for length estimation
            }
        }

        // Determine current column
        const colKey = `col${temp_col_state}_label_max`;

        // Measure the label length and update the max length for that column's category
        const label_length = label.length;
        if (label_max_dict[colKey][category] < label_length) {
            label_max_dict[colKey][category] = label_length;
        }
    });

    // Calculate the actual pixel widths for each column and category
    const categories = ['Class', 'LigandType', 'ReceptorFamily', 'Receptor'];
    const cols = ['Col1', 'Col2', 'Col3', 'Col4'];

    for (let i = 0; i < columns; i++) {
        const key = `col${i + 1}_label_max`;
        const max_label_values = label_max_dict[key];

        categories.forEach(category => {
            const columnKey = cols[i];

            if (max_label_values[category] !== -Infinity) {
                // Create a dummy text element to measure the width based on the max label length
                const dummyText = d3.select("body")
                    .append("svg")
                    .attr("class", "dummy-text")
                    .append("text")
                    .attr("font-size", styling_option[category].Fontsize) // Use the category as the key for styling_option
                    .attr("font-weight", styling_option[category].Bold ? "bold" : "normal")
                    .text("X".repeat(max_label_values[category]));

                const bbox = dummyText.node().getBBox();
                let estimatedLength = bbox.width * 0.8 + 20;

                if (category === 'Receptor') {
                    estimatedLength += 80; // Additional margin for Receptors
                }

                // Remove the dummy text element
                d3.select(".dummy-text").remove();

                // Update the spacing dict if this label is the largest so far for this column
                if (estimatedLength > Col_spacing_dict[columnKey]) {
                    Col_spacing_dict[columnKey] = estimatedLength;
                }
            }
        });
    }
    return Col_spacing_dict;
}

// Genreate the printing of the labels
function RenderListPlot_Labels(data, category_data, location, styling_option, Layout_dict, label_conversion_dict, label_names) {
    // ######################
    // ## Initialization   ##
    // ######################

    // Extract necessary values from the layout
    let columns = Layout_dict.columns;
    let col_max_label = Layout_dict.col_max_label;

    // Calculate total count
    const total_count = data.length;

    // Determine Col_break_number (automatic or custom)
    let Col_break_number = col_max_label === "Auto" ? Math.ceil(total_count / columns) : Layout_dict.Col_break_number;

    // ##################
    // ## Dimensions   ##
    // ##################
    let width = 0;
    const initial_height = 200;  // Starting height of the SVG
    const margin = { top: 40, right: 20, bottom: 20, left: 20 };

    // ## Calculate all spacing and dimensions ##
    let spacing_dict = Calculate_dimension(data, category_data, Col_break_number, columns, label_conversion_dict, label_names, styling_option);

    // Calculate total width based on spacing_dict and columns
    const col_list = ['Col1', 'Col2', 'Col3', 'Col4'];
    for (let i = 0; i < columns; i++) {
        width += spacing_dict[col_list[i]];
    }
    width += 45; // Add some padding

    // // Recalculate width if its less than legend bars
    let number_of_bar_legends = 0;
    Object.keys(Data_styling).forEach(function(column) {
        if (Data_styling[column].Data === "Yes" && Data_styling[column].Datatype === 'Continuous') {
            number_of_bar_legends++;
        }
    });
    let minimum_width = 0;
    if (number_of_bar_legends > 0) {
        minimum_width = 230*number_of_bar_legends+margin.left+margin.right;
        if (width < minimum_width) {
            width = minimum_width;
        }
    }

    // Initialize SVG container
    const svg = d3.select(`#${location}`)
        .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", initial_height + margin.top + margin.bottom)  // Start with initial height
        .attr("id", "ListPlot_visualization");

    // Set initial X & Y positions
    let yOffset = margin.top + 5 + 80;
    let yOffset_max = yOffset; // Track the maximum yOffset
    let xOffset = 0;

    // ######################
    // ## Label Handling   ##
    // ######################

    // Function to handle text rendering with special cases for IUPHAR and UniProt (only for Receptor category)
    function add_text(label, yOffset, category, label_key, bold = false, italic = false, underline = false, fontSize = '16px', color = 'black', label_names) {
        // Handle special cases for Receptors
        if (category === 'Receptor') {
            if (label_names === 'Protein') {

                const htmlEntities = {
                    '&alpha;': 'α', '&beta;': 'β', '&gamma;': 'γ', '&delta;': 'δ',
                    '&epsilon;': 'ε', '&zeta;': 'ζ', '&eta;': 'η', '&theta;': 'θ',
                    '&iota;': 'ι', '&kappa;': 'κ', '&lambda;': 'λ', '&mu;': 'μ',
                    '&nu;': 'ν', '&xi;': 'ξ', '&omicron;': 'ο', '&pi;': 'π',
                    '&rho;': 'ρ', '&sigma;': 'σ', '&tau;': 'τ', '&upsilon;': 'υ',
                    '&phi;': 'φ', '&chi;': 'χ', '&psi;': 'ψ', '&omega;': 'ω'
                };

                // Replace HTML entities with corresponding Unicode characters
                for (const [entity, char] of Object.entries(htmlEntities)) {
                    label = label.replace(new RegExp(entity, 'g'), char);
                }

                // Remove <i> tags and (-adrenoceptor| receptor) from the label
                label = label.replace(/<\/?i>/g, '')  // Remove <i> tags
                .replace(/(-adrenoceptor| receptor)/g, '');  // Remove adrenoceptor/receptor text


                // Create text element for IUPHAR receptor
                const textElement = svg.append('text')
                    .attr('x', margin.left + xOffset + 80)
                    .attr('y', yOffset)
                    .attr('class', category)
                    .style('dominant-baseline', 'middle') // Set vertical alignment
                    .style('font-weight', bold ? 'bold' : 'normal')
                    .style('font-style', italic ? 'italic' : 'normal')
                    .style('text-decoration', underline ? 'underline' : 'none')
                    .style('font-size', fontSize)
                    .style('fill', color);

                // Subscript handling for IUPHAR receptors
                const mainFontSize = parseFloat(fontSize);
                const subFontSize = mainFontSize * 0.75 + 'px';
                const parts = label.split(/(<sub>|<\/sub>)/);
                let inSub = false;

                // Handle subscripts
                parts.forEach(part => {
                    if (part === '<sub>') {
                        inSub = true;
                    } else if (part === '</sub>') {
                        inSub = false;
                    } else {
                        const tspan = textElement.append('tspan')
                            .text(part);

                        if (inSub) {
                            tspan.attr('dy', '0.3em') // Subscript positioning
                                .attr('font-size', subFontSize);
                        } else {
                            tspan.attr('dy', '-0.3em') // Reset positioning
                                .attr('font-size', fontSize);
                        }
                    }
                });

            } else if (label_names === 'Gene') {
                // Handle UniProt receptor labels
                label = label_conversion_dict[label_key]?.replace(/_human/g, '').toUpperCase() || label_key;
                svg.append('text')
                    .attr('x', margin.left + xOffset + 80)
                    .attr('y', yOffset)
                    .attr('class', category)
                    .attr('dy', '-0.3em') // Adjust this value to move the text higher
                    .style('dominant-baseline', 'middle') // Set vertical alignment
                    .style('font-weight', bold ? 'bold' : 'normal')
                    .style('font-style', italic ? 'italic' : 'normal')
                    .style('text-decoration', underline ? 'underline' : 'none')
                    .style('font-size', fontSize)
                    .style('fill', color)
                    .text(label);
            }
        } else if (category === 'ReceptorFamily') {
            // Handle ReceptorFamily (subscript handling but no label_names logic)
            label = label_key.replace(/( receptors|neuropeptide )/g, '').split(" (")[0];

            // Create text element for ReceptorFamily
            const textElement = svg.append('text')
                .attr('x', margin.left + xOffset)
                .attr('y', yOffset)
                .attr('class', category)
                .style('dominant-baseline', 'middle') // Set vertical alignment
                .style('font-weight', bold ? 'bold' : 'normal')
                .style('font-style', italic ? 'italic' : 'normal')
                .style('text-decoration', underline ? 'underline' : 'none')
                .style('font-size', fontSize)
                .style('fill', color);

            // Subscript handling for ReceptorFamily
            const mainFontSize = parseFloat(fontSize);
            const subFontSize = mainFontSize * 0.75 + 'px';
            const parts = label.split(/(<sub>|<\/sub>)/);
            let inSub = false;

            // Handle subscripts
            parts.forEach(part => {
                if (part === '<sub>') {
                    inSub = true;
                } else if (part === '</sub>') {
                    inSub = false;
                } else {
                    const tspan = textElement.append('tspan')
                        .text(part);

                    if (inSub) {
                        tspan.attr('dy', '0.3em') // Subscript positioning
                            .attr('font-size', subFontSize);
                    } else {
                        tspan.attr('dy', '-0.3em') // Reset positioning
                            .attr('font-size', fontSize);
                    }
                }
            });
        } else {
            // General label handling for other categories
            if (category === 'Class') {
                label = label_key.split(" (")[0]; // Trim Class
            }

            // Render text for non-receptor categories
            svg.append('text')
                .attr('x', margin.left + xOffset)
                .attr('y', yOffset)
                .attr('class', category)
                .style('dominant-baseline', 'middle') // Set vertical alignment
                .style('font-weight', bold ? 'bold' : 'normal')
                .style('font-style', italic ? 'italic' : 'normal')
                .style('text-decoration', underline ? 'underline' : 'none')
                .style('font-size', fontSize)
                .style('fill', color)
                .text(label);
        }
    }


    // ###########################
    // ## Render List of Labels ##
    // ###########################

    let label_counter = 1;
    let current_col = 1; // Track the current column being rendered

    for (let i = 0; i < data.length; i++) {
        const label = data[i];
        const category = category_data[i];
        const label_key = label; // Use original label before conversions
        const style = styling_option[category]; // Get the style options for the category

        // Add text for the label
        add_text(label, yOffset, category, label_key, style.Bold, style.Italic, style.Underline, style.Fontsize, style.Color,label_names);

        yOffset += 30; // Adjust Y-offset for the next label

        // Check if this is the new maximum yOffset encountered
        if (yOffset > yOffset_max) {
            yOffset_max = yOffset; // Update yOffset_max globally
        }

        label_counter++;

        // Handle column breaks when Col_break_number is reached
        if (label_counter > Col_break_number && current_col < columns) {
            current_col++;
            xOffset += spacing_dict[`Col${current_col - 1}`]; // Move to the next column
            label_counter = 1; // Reset label counter for the new column
            yOffset = margin.top + 5 + 80; // Reset yOffset for the new column
        }

        // Ensure the global maximum yOffset is tracked
        if (yOffset > yOffset_max) {
            yOffset_max = yOffset; // Update the maximum yOffset encountered
        }
    }

    // Rerender height of plot based on the max yOffset encountered
    svg.attr("height", yOffset_max + margin.top + margin.bottom);

    return svg;
}

// Handles the data visualization
function data_visualization(data, category_data, location, Layout_dict, data_styling, spacing_dict) {

    // ###########################
    // ## Initialize Variables  ##
    // ###########################

    // Set default margins, xOffset, yOffset, and columns based on Layout_dict
    const margin = { top: 40, right: 20, bottom: 20, left: 20 };
    let yOffset = margin.top + 5 + 80 + 5;
    let xOffset = 5;
    let columns = Layout_dict.columns;

    // Get the SVG container
    const svg = d3.select(`#${location} svg`);

    // Track the current column and Y-offsets
    let current_col = 1;
    let label_counter = 1;
    const Col_break_number = Layout_dict.Col_break_number || Math.ceil(data.length / columns);

    // Function to add different shapes
    function addShape(shapeType, x, y, size, fillColor) {
        switch (shapeType) {
            case 'circle':
                svg.append('circle')
                    .attr('cx', x)
                    .attr('cy', y)
                    .attr('r', size)
                    .style('dominant-baseline', 'middle') // Set vertical alignment
                    .style('stroke', 'black')
                    .style('stroke-width', 1)
                    .style('fill', fillColor);
                break;
            case 'rect':
                svg.append('rect')
                    .attr('x', x - size)
                    .attr('y', y - size)
                    .attr('width', size * 2)
                    .attr('height', size * 2)
                    .style('dominant-baseline', 'middle') // Set vertical alignment
                    .style('stroke', 'black')
                    .style('stroke-width', 1)
                    .style('fill', fillColor);
                break;
            case 'triangle':
                svg.append('path')
                    .attr('d', `M ${x} ${y - size} L ${x - size} ${y + size} L ${x + size} ${y + size} Z`)
                    .style('stroke', 'black')
                    .style('dominant-baseline', 'middle') // Set vertical alignment
                    .style('stroke-width', 1)
                    .style('fill', fillColor);
                break;
            case 'diamond':
                svg.append('path')
                    .attr('d', `M ${x} ${y - size} L ${x - size} ${y} L ${x} ${y + size} L ${x + size} ${y} Z`)
                    .style('stroke', 'black')
                    .style('dominant-baseline', 'middle') // Set vertical alignment
                    .style('stroke-width', 1)
                    .style('fill', fillColor);
                break;
            case 'star':
                const points = 5;
                const outerRadius = size;
                const innerRadius = size / 2.5;
                let starPath = '';
                for (let i = 0; i < points * 2; i++) {
                    const angle = (Math.PI / points) * i;
                    const radius = i % 2 === 0 ? outerRadius : innerRadius;
                    const xPoint = x + Math.cos(angle) * radius;
                    const yPoint = y + Math.sin(angle) * radius;
                    starPath += i === 0 ? `M ${xPoint} ${yPoint}` : `L ${xPoint} ${yPoint}`;
                }
                starPath += 'Z';
                svg.append('path')
                    .attr('d', starPath)
                    .style('stroke', 'black')
                    .style('dominant-baseline', 'middle') // Set vertical alignment
                    .style('stroke-width', 1)
                    .style('fill', fillColor);
                break;
            default:
                console.log('Unknown shape type');
        }
    }

    // Function to handle column breaks
    function handleColumnBreak() {
        if (label_counter > Col_break_number && current_col < columns) {
            current_col++;
            xOffset += spacing_dict[`Col${current_col - 1}`]; // Move to the next column
            label_counter = 1; // Reset label counter for the new column
            yOffset = margin.top + 5 + 80 + 5; // Reset yOffset for the new column
        }
    }

    // ###############################
    // ## Color Logic Function      ##
    // ###############################

    function getShapeColor(column, data_value) {
        const column_styling = data_styling[column];
        let color = 'black';

        if (column_styling.Datatype === 'Discrete') {
            const Color_list = ['black', 'red', 'blue', 'green'];
            if (data_value) {
                const valueString = String(data_value); // Ensures data_value is converted to a string
                color = Color_list.includes(valueString.toLowerCase()) ? valueString : 'black';
            } else {
                color = 'black';
            }
        } else if (column_styling.Datatype === 'Continuous') {
            const gradientScale = d3.scale.linear()
                .domain([column_styling.Data_min, column_styling.Data_max]);

            if (column_styling.data_color_complexity === 'One') {
                gradientScale.range(['#FFFFFF', column_styling.Data_color2]);
            } else if (column_styling.data_color_complexity === 'Two') {
                gradientScale.range([column_styling.Data_color1, column_styling.Data_color2]);
            } else if (column_styling.data_color_complexity === 'Three') {
                gradientScale.range([column_styling.Data_color1, '#FFFFFF', column_styling.Data_color2])
                    .domain([column_styling.Data_min, (column_styling.Data_min + column_styling.Data_max) / 2, column_styling.Data_max]);
            }

            color = gradientScale(data_value);
        }

        return color;
    }

    // ###############################
    // ## Iterate through data rows ##
    // ###############################

    for (let i = 0; i < data.length; i++) {
        const label = data[i];
        const category = category_data[i];

        // Check if category is 'Receptor' (only for receptors do we render shapes)
        if (category === 'Receptor') {
            // Check if data exists for the receptor (in the shape data)
            if (list_data_wow.hasOwnProperty(label)) {
                const receptorData = list_data_wow[label];

                // Col data checker
                const Col1_data_checker = data_styling.Col1.Data == "Yes";
                const Col2_data_checker = data_styling.Col2.Data == "Yes";
                const Col3_data_checker = data_styling.Col3.Data == "Yes";
                const Col4_data_checker = data_styling.Col4.Data == "Yes";

                // Col labels and shapes
                const Col1_shape = receptorData.hasOwnProperty('Value1') ? receptorData.Value1.toLowerCase() : false;
                const Col2_shape = receptorData.hasOwnProperty('Value3') ? receptorData.Value3.toLowerCase() : false;
                const Col3_shape = receptorData.hasOwnProperty('Value5') ? receptorData.Value5.toLowerCase() : false;
                const Col4_shape = receptorData.hasOwnProperty('Value7') ? receptorData.Value7.toLowerCase() : false;

                // Col data
                const Col1_data = receptorData.hasOwnProperty('Value2') ? receptorData.Value2 : false;
                const Col2_data = receptorData.hasOwnProperty('Value4') ? receptorData.Value4 : false;
                const Col3_data = receptorData.hasOwnProperty('Value6') ? receptorData.Value6 : false;
                const Col4_data = receptorData.hasOwnProperty('Value8') ? receptorData.Value8 : false;

                // Shapes and data rendering for each column
                const col1_XoffSet = 0;
                const col2_XoffSet = 20;
                const col3_XoffSet = 40;
                const col4_XoffSet = 60;

                const Shape_list = ['circle', 'rect', 'triangle', 'star', 'diamond'];

                // ### Column 1 ###
                if (Col1_data_checker && (Col1_shape || Col1_data)) {
                    const shape_color = Col1_data ? getShapeColor('Col1', Col1_data) : 'black';
                    addShape(Shape_list.includes(Col1_shape) ? Col1_shape : 'circle', margin.left + xOffset + col1_XoffSet, yOffset - 10, 6, shape_color);
                }

                // ### Column 2 ###
                if (Col2_data_checker && (Col2_shape || Col2_data)) {
                    const shape_color = Col2_data ? getShapeColor('Col2', Col2_data) : 'black';
                    addShape(Shape_list.includes(Col2_shape) ? Col2_shape : 'circle', margin.left + xOffset + col2_XoffSet, yOffset - 10, 6, shape_color);
                }

                // ### Column 3 ###
                if (Col3_data_checker && (Col3_shape || Col3_data)) {
                    const shape_color = Col3_data ? getShapeColor('Col3', Col3_data) : 'black';
                    addShape(Shape_list.includes(Col3_shape) ? Col3_shape : 'circle', margin.left + xOffset + col3_XoffSet, yOffset - 10, 6, shape_color);
                }

                // ### Column 4 ###
                if (Col4_data_checker && (Col4_shape || Col4_data)) {
                    const shape_color = Col4_data ? getShapeColor('Col4', Col4_data) : 'black';
                    addShape(Shape_list.includes(Col4_shape) ? Col4_shape : 'circle', margin.left + xOffset + col4_XoffSet, yOffset - 10, 6, shape_color);
                }
            }

            // Increment yOffset for the next label and shape
            yOffset += 30;
            label_counter++;

            // Handle column break
            handleColumnBreak();

        } else {
            // If the category is not 'Receptor', simply move the Y offset without rendering shapes
            yOffset += 30;
            label_counter++;

            // Handle column break
            handleColumnBreak();
        }
    }

    // #############################
    // ## Render Legend Bars on Top ##
    // #############################

    let bar_index = 0;
    Object.keys(data_styling).forEach(function(column) {
        if (data_styling[column].Data === "Yes" && data_styling[column].Datatype === 'Continuous') {
            const legendWidth = 200; // Width of the legend bar
            const data_fontsize = 14; // Adjust as needed
            const lowest_value = data_styling[column].Data_min;
            const highest_value = data_styling[column].Data_max;
            const spacing_bar = 30;
            const bar_height = 20;
            const text_off_set = 35;
            const x_off_set = 100;
            const y_off_set = 20;

            // Calculate the x position for the current bar
            const x_position = bar_index * (legendWidth + spacing_bar) + x_off_set;

            // Increment bar index for the next bar
            bar_index++;

            const legend_svg = svg.append("g"); // Append group for the legend bar
            const gradientId = `Gradient_${column}`;
            const defs = svg.append('defs');

            // Add centered text on top of the bar specifying the data column
            legend_svg.append("text")
                .attr('x', x_position + legendWidth / 2) // Center the text horizontally over the bar
                .attr('y', bar_height - 10 + y_off_set) // Adjust Y position as needed, e.g., 10px above the bar
                .style("font-size", `${data_fontsize}px`)
                .style("font-family", "sans-serif")
                .style("text-anchor", "middle") // Center the text horizontally
                .text(`Data column ${column.substr(3)}`); // Use the column name as the text

            // Gradient definition
            const gradient = defs.append('linearGradient')
                .attr('id', gradientId)
                .attr('x1', '0%')
                .attr('x2', '100%')
                .attr('y1', '0%')
                .attr('y2', '0%');

            // Adjust gradient stops based on color complexity
            if (data_styling[column].data_color_complexity === 'One') {
                // One color: white to Data_color1
                gradient.append('stop')
                    .attr('offset', '0%')
                    .attr('stop-color', '#FFFFFF'); // Start with white

                gradient.append('stop')
                    .attr('offset', '100%')
                    .attr('stop-color', data_styling[column].Data_color2); // End with color1
            } else if (data_styling[column].data_color_complexity === 'Two') {
                // Two colors: Data_color1 to Data_color2
                gradient.append('stop')
                    .attr('offset', '0%')
                    .attr('stop-color', data_styling[column].Data_color1); // Start with color1

                gradient.append('stop')
                    .attr('offset', '100%')
                    .attr('stop-color', data_styling[column].Data_color2); // End with color2
            } else if (data_styling[column].data_color_complexity === 'Three') {
                // Three colors: Data_color1 to white in the middle, then to Data_color2
                gradient.append('stop')
                    .attr('offset', '0%')
                    .attr('stop-color', data_styling[column].Data_color1); // Start with color1

                gradient.append('stop')
                    .attr('offset', '50%')
                    .attr('stop-color', '#FFFFFF'); // Middle with white

                gradient.append('stop')
                    .attr('offset', '100%')
                    .attr('stop-color', data_styling[column].Data_color2); // End with color2
            }

            // Border around the rectangle
            legend_svg.append('rect')
                .attr('x', x_position) // Align with the gradient rectangle
                .attr('y', bar_height + y_off_set) // Align with the gradient rectangle
                .attr('width', legendWidth) // Same width as the gradient rectangle
                .attr('height', bar_height) // Same height as the gradient rectangle
                .style('fill', 'none')
                .style('stroke', 'black')
                .style('stroke-width', 2);

            // Rectangle for the gradient
            legend_svg.append('rect')
                .attr('x', x_position) // Horizontal spacing between bars
                .attr('y', bar_height + y_off_set) // Adjust Y position as needed
                .attr('width', legendWidth) // Legend bar width
                .attr('height', bar_height) // Legend bar height
                .style('fill', `url(#${gradientId})`);

            // Minimum value text
            legend_svg.append("text")
                .attr('x', x_position) // Align with the bar
                .attr('y', bar_height + text_off_set + y_off_set)
                .style("font-size", `${data_fontsize}px`)
                .style("font-family", "sans-serif")
                .text(lowest_value);

            // Middle value text (if applicable)
            if (data_styling[column].data_color_complexity === 'Three') {
                legend_svg.append("text")
                    .attr('x', x_position + legendWidth / 2)
                    .attr('y', bar_height + text_off_set + y_off_set)
                    .style("font-size", `${data_fontsize}px`)
                    .style("font-family", "sans-serif")
                    .style("text-anchor", "middle")
                    .text((highest_value + lowest_value) / 2); // Middle value
            }

            // Maximum value text
            legend_svg.append("text")
                .attr('x', x_position + legendWidth)
                .attr('y', bar_height + text_off_set + y_off_set)
                .style("font-size", `${data_fontsize}px`)
                .style("font-family", "sans-serif")
                .style("text-anchor", "end")
                .text(highest_value);
        }
    });
}

// #################
// ###  HEATMAP  ###
// #################

// Initialize the datastyling dict
function heatmap_DataStyling() {
    var heatmap_DataStyling = {Number_of_colors: "Three",
        min_color: "#1a80bb",
        middle_color: "#FFFFFF",
        max_color: "#a00000",
        rotation: 90,
        label_position: 'Bottom',
        label_fontsize: 12,
        receptor_fontsize: 12,
        datalabels: true,
        data_border: true,
        data_fontsize: 12,
        legend_label: "Value intensity",
        LabelType: 'UniProt'

    }
    return heatmap_DataStyling
}

// Function to handle row labels based on label type
function handleRowLabels(textElement, label, labelType, fontSize) {
    const htmlEntities = {
        '&alpha;': 'α',
        '&beta;': 'β',
        '&gamma;': 'γ',
        '&delta;': 'δ',
        '&epsilon;': 'ε',
        '&zeta;': 'ζ',
        '&eta;': 'η',
        '&theta;': 'θ',
        '&iota;': 'ι',
        '&kappa;': 'κ',
        '&lambda;': 'λ',
        '&mu;': 'μ',
        '&nu;': 'ν',
        '&xi;': 'ξ',
        '&omicron;': 'ο',
        '&pi;': 'π',
        '&rho;': 'ρ',
        '&sigma;': 'σ',
        '&tau;': 'τ',
        '&upsilon;': 'υ',
        '&phi;': 'φ',
        '&chi;': 'χ',
        '&psi;': 'ψ',
        '&omega;': 'ω'
    };

    // Handle label transformations based on labelType
    let transformedLabel = label; // Initialize transformedLabel with the original label

    if (labelType === 'UniProt') {
        transformedLabel = label.replace(/_human/g, '').toUpperCase();
        textElement.text(transformedLabel);
    } else if (labelType === 'IUPHAR') {

        // Retrieve the IUPHAR label from your converter object
        transformedLabel = label_converter.UniProt_to_IUPHAR_converter[label];
        // Replace HTML entities with corresponding Unicode characters
        for (const [entity, char] of Object.entries(htmlEntities)) {
            transformedLabel = transformedLabel.replace(new RegExp(entity, 'g'), char);
        }
        // Remove <i> and </i> tags from the label
        transformedLabel = transformedLabel.replace(/<\/?i>/g, '');
        // remove receptor and adrenoceptor
        transformedLabel = transformedLabel.replace("-adrenoceptor", '').replace(" receptor", '');
        // Clear existing content in textElement if needed
        textElement.text('')

        // Additional transformations specific to IUPHAR labels
        const parts = transformedLabel.split(/(<sub>|<\/sub>)/); // Split label into parts including <sub> tags
        let inSub = false; // Flag to track if we are inside a subscript

        // Calculate the subscript font size (e.g., 75% of the main font size)
        const mainFontSize = parseFloat(fontSize);
        const subFontSize = parseInt(mainFontSize * 0.75,10);

        // Handle each part of the label
        parts.forEach(part => {
            if (part === '<sub>') {
                inSub = true;
            } else if (part === '</sub>') {
                inSub = false;
            } else {
                const tspan = textElement.append('tspan').text(part).style("font-family", "sans-serif");

                // Adjust dy based on subscript flag
                if (inSub) {
                    tspan.attr('dy', '0.5em').attr('font-size', subFontSize); // Adjust dy for subscript
                } else {
                    tspan.attr('dy', '0.1em').attr('font-size', mainFontSize); // Default dy for regular text
                }
            }
        });
    }
}

// Create the heatmap
function Heatmap(data, location, heatmap_DataStyling,label_x_converter) {

    const margin = { top: 30, right: 100, bottom: 30, left: 60 }; // Adjusted margin for row labels
    const rows = Object.keys(data);
    const cols = Object.keys(data[rows[0]]);
    const col_labels = cols.map(col => label_x_converter[col]);

    const rotation = heatmap_DataStyling.rotation;
    const data_labels = heatmap_DataStyling.datalabels;
    const labelType = heatmap_DataStyling.LabelType;
    const label_position = heatmap_DataStyling.label_position;
    const label_fontsize = heatmap_DataStyling.label_fontsize;
    const receptor_fontsize = heatmap_DataStyling.receptor_fontsize;
    const data_fontsize = heatmap_DataStyling.data_fontsize;
    let legend_y_position = 30;

    let rowLabelWidth;
    const baseWidth = 20;  // Minimum column width
    const chartData = [];
    let highest_value = -Infinity;
    let lowest_value = Infinity;

    rows.forEach(row => {
      cols.forEach(col => {
        const value = data[row][col];
        chartData.push({ row, col, value });
        if (value > highest_value) {
          highest_value = value;
        }
        if (value < lowest_value) {
            lowest_value = value;
        }
      });
    });

    // Calculate the length of the longest data value

    const longestDataValue = d3.max(chartData, d => Number(parseFloat(d.value).toFixed(1)).toString().length);

    // Calculate width based on the longest data value
    const calculatedDataValueWidth = longestDataValue * 40 * (label_fontsize / 14);

    if (rotation === 90 || rotation === 45) {
        // Use the larger value between base width and calculated data value width for rotation
        rowLabelWidth = Math.max(baseWidth, calculatedDataValueWidth);
    } else {
        // Calculate the width based on the longest column label
        const calculatedLabelWidth = d3.max(col_labels, d => d.length * 40 * (label_fontsize / 14));

        // Choose the maximum value between base width, calculated label width, and calculated data value width
        rowLabelWidth = Math.max(baseWidth, calculatedLabelWidth, calculatedDataValueWidth);
    }

    // Calculate the longest column label length
    const longestLabel = d3.max(col_labels, d => d.length);
    if (label_position === 'Top') {
        // Calculate the longest column label length from the data
        const fontSize = label_fontsize; // Assuming a font size of 14px
        // Adjust margin calculation with a scaling factor
        if (rotation === 90) {
            margin.top = longestLabel * (fontSize * 0.6); // Adding extra space for padding
        }
        legend_y_position = 0;
    } else if (label_position === 'Bottom') {
        if (rotation === 90) {
             // Calculate the longest column label length from the data
            const fontSize = label_fontsize; // Assuming a font size of 14px
            legend_y_position = longestLabel * (fontSize * 0.55)
        }
    }

    const width = (20 * cols.length) + rowLabelWidth;
    const height = (20 * rows.length);
    const svg_home = d3.select("#" + location)
      .append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + (rotation === 90 ? (margin.top * longestLabel/3) : (margin.top*2))) // Needs to account for label length or something like it.
      .attr("id", "Heatmap_plot_svg");

    //  x / y scale transformers
    const x = d3.scale.ordinal()
        .rangeBands([0, width])
        .domain(cols);

    const y = d3.scale.ordinal()
      .rangeBands([height, 0])
      .domain(rows.reverse());

    if (heatmap_DataStyling.Number_of_colors === 'Three') {
        var myColor = d3.scale.linear()
        .range([heatmap_DataStyling.min_color, heatmap_DataStyling.middle_color, heatmap_DataStyling.max_color]) // Adjust colors here
        .domain([lowest_value, (highest_value + lowest_value) / 2, highest_value]); // Adjust domain based on your data
    } else if (heatmap_DataStyling.Number_of_colors === 'One' || heatmap_DataStyling.Number_of_colors === 'Two') {
        var myColor = d3.scale.linear()
        .range([heatmap_DataStyling.min_color, heatmap_DataStyling.max_color]) // Adjust colors here
        .domain([lowest_value, highest_value]); // Adjust domain based on your data
    }

    const svg = svg_home.append("g")
        .attr("transform", `translate(${margin.left * 2}, ${margin.top})`);

    let xAxis;
    if (label_position === 'Bottom') {
    xAxis = svg.append("g")
        .attr("transform", `translate(0, ${height})`)
        .attr('id', 'Xaxis')
        .call(d3.svg.axis().scale(x).orient("bottom").tickSize(0));
    } else if (label_position === 'Top') {
    xAxis = svg.append("g")
        .attr("transform", `translate(0, -30)`)
        .attr('id', 'Xaxis')
        .call(d3.svg.axis().scale(x).orient("top").tickSize(0));
    }

    xAxis.selectAll("text")
    .attr("class", "column-label")
    .style("font-size", `${label_fontsize}px`)
    .style("font-family", "sans-serif")
    .each(function(d, i) {
        const col = cols[i]; // Assuming cols is defined somewhere in your code
        const text = d3.select(this);

        // Set text content based on label_x_converter
        text.text(label_x_converter[col]);

        if (rotation === 0) {
            text.style("text-anchor", "middle")
                .attr("transform", "rotate(0)")
                .attr("dx", `0`)
                .attr("dy", `${label_fontsize * 1.5}px`); // Adjust dy based on font size
        } else if (rotation === 45) {
            text.style("text-anchor", label_position === 'Bottom' ? "end" : "start")
                .attr("transform", "rotate(-45)")
                .attr("dx", label_position === 'Bottom' ? `${-0.5 * label_fontsize}px` : `${-1.1 * label_fontsize}px`)
                .attr("dy", label_position === 'Bottom' ? `${0.5 * label_fontsize}px` : `${1.1 * label_fontsize}px`);
        } else if (rotation === 90) {
            text.style("text-anchor", label_position === 'Bottom' ? "end" : "start")
                .attr("transform", "rotate(-90)")
                .attr("dx", label_position === 'Bottom' ? `${-0.3 * label_fontsize}px` : `${-1.3 * label_fontsize}px`)
                .attr("dy", label_position === 'Bottom' ? `0px` : `${0.6 * label_fontsize}px`);
        }
    });

    // Y labels
    const yAxis = svg.append("g")
        .attr('id', 'Yaxis')
        .call(d3.svg.axis().scale(y).orient("left").tickSize(0))
        .selectAll("text")
        .style("font-size", `${receptor_fontsize}px`)
        .style("font-family", "sans-serif");

    yAxis.each(function(d) {
        const textElement = d3.select(this);

        // Call function to handle row labels based on label type
        handleRowLabels(textElement, d, labelType,receptor_fontsize);
    });

    // Heatmap data boxes

    const rects = svg.selectAll("rect")
      .data(chartData, d => `${d.row}:${d.col}`)
      .enter()
      .append("rect")
      .attr("x", d => x(d.col))
      .attr("y", d => y(d.row))
      .attr("width", x.rangeBand())
      .attr("height", y.rangeBand())
      .style("fill", d => myColor(d.value))
      .each(function(d) {
        if (heatmap_DataStyling.data_border) {
            d3.select(this)
                .style('stroke', 'black')
                .style('stroke-width', 0.5);
        }
    });


    if (data_labels) {
      rects.each(function(d) {
        const rect = d3.select(this);
        const textColor = getContrastColor(myColor(d.value));

        svg.append("text")
          .attr("x", +rect.attr("x") + rect.attr("width") / 2)
          .attr("y", +rect.attr("y") + rect.attr("height") / 2)
          .attr("dy", ".35em")
          .attr("text-anchor", "middle")
          .style("font-size", `${data_fontsize}px`)
          .style("font-family", "sans-serif")
          .style("fill", textColor)
          .text(Number(parseFloat(d.value).toFixed(1)));  // Round and fix to 1 decimal place
      });
    }

    function getContrastColor(hexColor) {
      const rgb = d3.rgb(hexColor);
      const brightness = (rgb.r * 0.299 + rgb.g * 0.587 + rgb.b * 0.114);
      return brightness > 128 ? "black" : "white";
    }

    // ##########
    // # Legend #
    // ##########

    // Adjust for the width of the legend to center it
    let legendX;
    if (cols.length % 2 === 0) {
        // Even number of columns
        const leftMiddleColIndex = (cols.length / 2) - 1;
        const rightMiddleColIndex = cols.length / 2;
        const leftMiddleX = x(cols[leftMiddleColIndex]) + x.rangeBand() / 2;
        const rightMiddleX = x(cols[rightMiddleColIndex]) + x.rangeBand() / 2;
        legendX = (leftMiddleX + rightMiddleX) / 2;
    } else {
        // Odd number of columns
        const middleColIndex = Math.floor(cols.length / 2);
        legendX = x(cols[middleColIndex]) + x.rangeBand() / 2; // Center of the middle column
    }

    // Adjust for the width of the legend to center it
    let legendWidth;
    if (width <= 300) {
        legendWidth = width;
    } else {
        legendWidth = 300;
    }
    const adjustedLegendX = legendX - (legendWidth / 2);

    // Position the legend using the adjustedLegendX
    const legend_svg = svg.append("g")
        .attr("transform", `translate(${adjustedLegendX}, ${height + legend_y_position})`);

    const gradientId = "Heatmap_gradient";
    const defs = svg_home.append('defs');

    // Gradient definition
    const gradient = defs.append('linearGradient')
      .attr('id', gradientId)
      .attr('x1', '0%')
      .attr('x2', '100%')
      .attr('y1', '0%')
      .attr('y2', '0%');

    gradient.append('stop')
      .attr('offset', '0%')
      .attr('stop-color', heatmap_DataStyling.min_color); // Start with red at 0%

    if (heatmap_DataStyling.Number_of_colors === 'Three') {
        gradient.append('stop')
        .attr('offset', '50%')
        .attr('stop-color', heatmap_DataStyling.middle_color); // Middle point is white
    }
    gradient.append('stop')
      .attr('offset', '100%')
      .attr('stop-color', heatmap_DataStyling.max_color); // End with green at 100%

    const gradientRect = legend_svg.append('rect')
      .attr('x', 1) // Adjust for border
      .attr('y', 16) // Adjust for border
      .attr('width', legendWidth - 2) // Adjust for border
      .attr('height', 14) // Adjust for border
      .style('fill', `url(#${gradientId})`);

    if (heatmap_DataStyling.data_border) {
      const borderRect = legend_svg.append('rect')
      .attr('x', 0)
      .attr('y', 15)
      .attr('width', legendWidth)
      .attr('height', 15)
      .style('fill', 'none')
      .style('stroke', 'black')
      .style('stroke-width', 1);
    }
    legend_svg.append("text")
      .attr('x', 0)
      .attr('y', 50)
      .style("font-size", `${data_fontsize}px`)
      .style("font-family", "sans-serif")
      .text(Number(parseFloat(lowest_value).toFixed(1)));

    if (heatmap_DataStyling.Number_of_colors === 'Three') {
        legend_svg.append("text")
        .attr('x', legendWidth / 2)
        .attr('y', 50)
        .style("font-size", `${data_fontsize}px`)
        .style("font-family", "sans-serif")
        .style("text-anchor", "middle")
        .text(Number(parseFloat((highest_value + lowest_value) / 2).toFixed(1)));
    }

    legend_svg.append("text")
      .attr('x', legendWidth)
      .attr('y', 50)
      .style("font-size", `${data_fontsize}px`)
      .style("font-family", "sans-serif")
      .style("text-anchor", "end")
      .text(Number(parseFloat(highest_value).toFixed(1)));

    // Rerender height of plot as the last thing
    let label_length_final;
    if (rotation === 90) {
        label_length_final = Math.ceil(-128.69+6.61*longestLabel+11.19*label_fontsize)
        svg_home.attr("height",height + margin.bottom + label_length_final+55);
    } else {
        svg_home.attr("height",height + margin.bottom + 55 + 55);
    }

  }

// #################
// ###  GPCRome  ###
// #################

function GPCRome_initializeOdorantData(data) {
    // Initialize the GPCRomes
    let GPCRomes = {
        GPCRome_O1: {},
        GPCRome_O2_ext: {},
        GPCRome_O2_mid: {},
        GPCRome_O2_int: {},
    };

    // Helper function to add items to the GPCRome using receptor family as key
    function addItemsToCircle(GPCRome, items) {
        Object.keys(items).forEach(ligandType => {
            const receptorFamilies = items[ligandType];

            // Debugging - log the receptorFamilies structure
            // console.log(`Processing ligandType: ${ligandType}`);
            // console.log(`Receptor families:`, receptorFamilies);
            // Ensure receptorFamilies is an object
            if (typeof receptorFamilies === 'object' && receptorFamilies !== null) {
                Object.keys(receptorFamilies).forEach(family => {
                    if (!GPCRome[family]) {
                        GPCRome[family] = [];
                    }

                    const receptors = receptorFamilies[family];
                    if (Array.isArray(receptors)) {
                        GPCRome[family].push(...receptors);
                    } else {
                        console.warn(`Unexpected data format for family ${family}`);
                    }
                });
            } else {
                console.warn(`Unexpected data format for ligand type ${ligandType}`);
            }
        });
    }

    // Process the data
    Object.keys(data).forEach(className => {
        const classData = data[className];
        // Handle cases where classData is not an object
        if (typeof classData !== 'object' || classData === null) {
            console.warn(`Expected object but got ${typeof classData} for class ${className}`);
            return;
        }

        // Circle 1 : "Class O2 (tetrapod specific odorant) EXT"
        if (className === "Class O2 (tetrapod specific odorant) EXT") {
            addItemsToCircle(GPCRomes.GPCRome_O2_ext, classData);
        }

        // Circle 2 : "Class O2 (tetrapod specific odorant) MID"
        if (className === "Class O2 (tetrapod specific odorant) MID") {
            addItemsToCircle(GPCRomes.GPCRome_O2_mid, classData);
        }

        // Circle 3 : "Class O2 (tetrapod specific odorant) INT"
        if (className === "Class O2 (tetrapod specific odorant) INT") {
            addItemsToCircle(GPCRomes.GPCRome_O2_int, classData);
        }

        // Circle 4: "Class O1 (fish-like odorant)"
        if (className === "Class O1 (fish-like odorant)") {
            addItemsToCircle(GPCRomes.GPCRome_O1, classData);
        }

    });

    // Convert the arrays to unique values
    Object.keys(GPCRomes).forEach(GPCRomeKey => {
        Object.keys(GPCRomes[GPCRomeKey]).forEach(familyKey => {
            GPCRomes[GPCRomeKey][familyKey] = Array.from(new Set(GPCRomes[GPCRomeKey][familyKey]));
        });
    });

    return GPCRomes;
}

// Initialize the data for the handle of the GPCRome from the view
function GPCRome_initializeData(data) {
    // Initialize the GPCRomes
    let GPCRomes = {
        GPCRome_A: {},
        GPCRome_AO: {},
        GPCRome_B1: {},
        GPCRome_B2: {},
        GPCRome_C: {},
        GPCRome_F: {},
        GPCRome_T: {},
        GPCRome_CL: {}
    };

    // Helper function to add items to the GPCRome using receptor family as key
    function addItemsToCircle(GPCRome, items) {
        Object.keys(items).forEach(ligandType => {
            const receptorFamilies = items[ligandType];

            // Ensure receptorFamilies is an object
            if (typeof receptorFamilies === 'object' && receptorFamilies !== null) {
                Object.keys(receptorFamilies).forEach(family => {
                    if (!GPCRome[family]) {
                        GPCRome[family] = [];
                    }

                    const receptors = receptorFamilies[family];
                    if (Array.isArray(receptors)) {
                        GPCRome[family].push(...receptors);
                    } else {
                        console.warn(`Unexpected data format for family ${family}`);
                    }
                });
            } else {
                console.warn(`Unexpected data format for ligand type ${ligandType}`);
            }
        });
    }

    // Process the data
    Object.keys(data).forEach(className => {
        const classData = data[className];

        // Handle cases where classData is not an object
        if (typeof classData !== 'object' || classData === null) {
            console.warn(`Expected object but got ${typeof classData} for class ${className}`);
            return;
        }

        // Circle 1: All of "Class A (Rhodopsin)" excluding "Orphan receptors"
        if (className === "Class A (Rhodopsin)") {
            Object.keys(classData).forEach(ligandType => {
                if (ligandType !== "Orphan receptors" && ligandType !== "Olfactory receptors") {
                    addItemsToCircle(GPCRomes.GPCRome_A, { [ligandType]: classData[ligandType] });
                }
            });
        }

        // Circle 2: "Orphan receptors" from "Class A (Rhodopsin)"
        if (className === "Class A (Rhodopsin)" && classData["Orphan receptors"]) {
            addItemsToCircle(GPCRomes.GPCRome_AO, { "Orphan receptors": classData["Orphan receptors"] });
        }

        // Circle 3: "Class B1 (Secretin)" or "Class B2 (Adhesion)"
        if (className === "Class B1 (Secretin)") {
            addItemsToCircle(GPCRomes.GPCRome_B1, classData);
        }

        // Circle 3: "Class B1 (Secretin)" or "Class B2 (Adhesion)"
        if (className === "Class B2 (Adhesion)") {
            addItemsToCircle(GPCRomes.GPCRome_B2, classData);
        }

        // Circle 4: "Class C (Glutamate)"
        if (className === "Class C (Glutamate)") {
            addItemsToCircle(GPCRomes.GPCRome_C, classData);
        }

        // Circle 5: "Class F (Frizzled)"
        if (className === "Class F (Frizzled)") {
            addItemsToCircle(GPCRomes.GPCRome_F, classData);
        }

        // Circle 6: "Class T (Taste 2)"
        if (className === "Class T2 (Taste 2)") {
            addItemsToCircle(GPCRomes.GPCRome_T, classData);
        }

        // Circle 7: "Other GPCRs"
        if (className === "Other GPCRs") {
            addItemsToCircle(GPCRomes.GPCRome_CL, classData);
        }
    });

    // Convert the arrays to unique values
    Object.keys(GPCRomes).forEach(GPCRomeKey => {
        Object.keys(GPCRomes[GPCRomeKey]).forEach(familyKey => {
            GPCRomes[GPCRomeKey][familyKey] = Array.from(new Set(GPCRomes[GPCRomeKey][familyKey]));
        });
    });

    return GPCRomes;
}

// Reformat the labels (manual curated)
function GPCRome_formatTextWithHTML(text, Family_list) {
    // Apply all the replacements step by step
    let formattedText = text
        .replace(" receptors", '')
        .replace(" receptor", '')
        .replace("-adrenoceptor", '')
        .replace(" receptor-", '-')
        .replace("<sub>", '</tspan><tspan baseline-shift="-20%">')
        .replace("</sub>", '</tspan><tspan>')
        .replace("<i>", '</tspan><tspan font-style="italic">')
        .replace("</i>", '</tspan><tspan>')
        .replace("Long-wave-sensitive", 'LWS')
        .replace("Medium-wave-sensitive", 'MWS')
        .replace("Short-wave-sensitive", 'SWS')
        .replace("Olfactory", 'OLF')
        .replace("calcitonin-like receptor", 'CLR')
        .replace("5-Hydroxytryptamine", '5-HT');

    // Capitalize the first letter after applying all replacements
    function capitalizeFirstLetter(str) {
        let match = str.match(/^[^a-zA-Z]*([a-zA-Z])/);
        if (match) {
            let index = match.index + match[0].length - 1;
            let firstLetter = str[index].toUpperCase(); // Capitalize the first letter
            return str.slice(0, index) + firstLetter + str.slice(index + 1);
        }
        return str; // If no alphabetic characters, return unchanged
    }

    // Apply capitalization after all replacements
    formattedText = capitalizeFirstLetter(formattedText);

    // Check if the text is in the Family_list for additional formatting
    const isInFamilyList = Family_list.includes(text);

    // Apply additional formatting if the text is in the Family_list
    if (isInFamilyList) {
        formattedText = formattedText
            .replace(/( receptors|neuropeptide )/g, '') // Remove specific substrings
            .replace(/(-releasing)/g, '-rel.') // Abbreviate specific substrings
            .replace(/(-concentrating)/g, '-conc.') // Abbreviate specific substrings
            .replace(/( and )/g, ' & ') // Replace "and" with "&"
            .replace(/(GPR18, GPR55 & GPR119)/g, 'GPR18, 55 & 119') // Special case formatting
            .replace(/(Class C Orphans)/g, 'Orphans') // Replace "Class C Orphans"
            .split("</tspan>")[0] // Keep only part before the first closing tspan tag
            .split(" (")[0]; // Keep only the part before the first " (" parenthesis
    }

    return formattedText;
}


// Draw / generate the GPCRome plot
function Draw_GPCRomes(layout_data, fill_data, location, GPCRome_styling, odorant = false) {

    dimensions = { height: 1000, width: 1000 };

    const Spacing = GPCRome_styling.Spacing;
    const datatype = GPCRome_styling.datatype;
    const family = GPCRome_styling.family
    const showIcon = GPCRome_styling.showIcon;  // Get the icon visibility state
    const Font_family = "Arial"
    const Class_fontsize = "20px"
    const Receptor_and_family_fontsize = "11px"

    const svg = d3v4.select("#" + location)
    .append("svg")
    .attr("id", location+'_svg')  // Add the id here for download reference
    .attr("width", dimensions.width)
    .attr("height", dimensions.height)
    .attr("xmlns", "http://www.w3.org/2000/svg")  // Add the SVG namespace
    .attr("xmlns:xlink", "http://www.w3.org/1999/xlink");  // Add the xlink namespace for images

   // Get the image URL from the data-attribute
    const imageUrl = document.getElementById("image-container").getAttribute("data-image-url");

    if (showIcon) {
        var img = new Image();
        img.crossOrigin = "Anonymous";  // Ensure cross-origin handling
        img.src = imageUrl;

        img.onload = function() {
            var canvas = document.createElement('canvas');
            canvas.width = img.width;
            canvas.height = img.height;
            var context = canvas.getContext('2d');
            context.drawImage(img, 0, 0);

            // Convert the image to a base64-encoded string
            var dataUrl = canvas.toDataURL('image/png');

            // Append the image to the SVG using 'xlink:href' (D3 v4 compatible)
            svg.append("image")
                .attr("xlink:href", dataUrl)  // Use 'xlink:href' for D3 v4 compatibility
                .attr("x", 0)  // Top-left corner
                .attr("y", -45)  // Top-left corner
                .attr("width", 230)  // Set width for the image
                .attr("height", 230)  // Set height for the image
                .attr("class", "toggle-image");  // Add a class to control visibility
        };
    }

    // SORT the data

    // Function to perform natural sorting
    const naturalSort = (obj) => {
        const collator = new Intl.Collator(undefined, { numeric: true, sensitivity: 'base' });

        // Get the keys of the object and sort them naturally
        const sortedKeys = Object.keys(obj).sort(collator.compare);

        // Create a new sorted object based on the sorted keys
        const sortedObj = {};
        sortedKeys.forEach(key => {
            sortedObj[key] = obj[key];  // Maintain the original data under sorted keys
        });

        return sortedObj;
    };

    // Iterate over each GPCRome_X object and sort its receptor families (keys)
    Object.keys(layout_data).forEach(dartKey => {
        layout_data[dartKey] = naturalSort(layout_data[dartKey]);
    });

    // Now call Draw_a_GPCRome for both updated GPCRome_A and GPCRome_AO

    // First Draw_a_GPCRome with GPCRome_A, now without the two receptors

    if (odorant) {

        Draw_a_GPCRome(layout_data.GPCRome_O2_ext, fill_data, 0, dimensions, Spacing, odorant);
        Draw_a_GPCRome(layout_data.GPCRome_O2_mid, fill_data, 1, dimensions, Spacing, odorant);
        Draw_a_GPCRome(layout_data.GPCRome_O2_int, fill_data, 2, dimensions, Spacing, odorant);
        Draw_a_GPCRome(layout_data.GPCRome_O1, fill_data, 3, dimensions, Spacing, odorant);

    } else {

      // Number of last entries to transfer and remove
      const N = 16;  // Change this value to 2, 3, or any number you want

      // Get the keys of the GPCRome_A object
      const dartAKeys = Object.keys(layout_data.GPCRome_A);

      // Get the last N keys
      const lastNKeys = dartAKeys.slice(-N);

      // Create a copy of GPCRome_AO and add the new entries from GPCRome_A
      const updatedGPCRome_AO = {
          ...lastNKeys.reduce((acc, key) => {
              acc[key] = layout_data.GPCRome_A[key]; // Add the last N entries from GPCRome_A
              return acc;
          }, {}),
          ...layout_data.GPCRome_AO // Spread the original GPCRome_AO entries
      };

      // Create a new GPCRome_A object that excludes the last N entries
      const updatedGPCRome_A = {
          ...layout_data.GPCRome_A
      };

      // Remove the last N properties from GPCRome_A
      lastNKeys.forEach(key => {
          delete updatedGPCRome_A[key];
      });

      Draw_a_GPCRome(updatedGPCRome_A, fill_data, 0, dimensions, Spacing);

      // Second Draw_a_GPCRome with GPCRome_AO, now with the two receptors added
      Draw_a_GPCRome(updatedGPCRome_AO, fill_data, 1, dimensions, Spacing);

      Draw_a_GPCRome({...layout_data.GPCRome_B1, ...layout_data.GPCRome_B2}, fill_data, 2, dimensions,Spacing);
      Draw_a_GPCRome({...layout_data.GPCRome_C, ...layout_data.GPCRome_F}, fill_data, 3, dimensions,Spacing);
      Draw_a_GPCRome({...layout_data.GPCRome_T,...layout_data.GPCRome_CL}, fill_data, 4, dimensions,Spacing);

    }
    function Draw_a_GPCRome(label_data, fill_data, level, dimensions, Spacing, odorant = false) {

        // Define SVG dimensions
        const width = dimensions.width;
        const height = dimensions.height;
        const label_offset = 7; // Increased offset to push labels outward
        let GPCRome_radius;

        if (odorant) {
          GPCRome_radius = Math.min(width, height) / 2 - 60 - ((level === 3) ? (100 * level) : (95 * level)); // Radius for each GPCRome
        } else {
          GPCRome_radius = Math.min(width, height) / 2 - 60 - ((level === 4) ? (90 * level) : (85 * level)); // Radius for each GPCRome
        }

        let values = [];
        let Family_list  = []
        let Family_exclude = ['Other GPCR orphans','Taste 2 receptors']
        if (family) {
            if (Spacing && Object.keys(label_data).length > 1) {
                // If Spacing is true and there are multiple keys
                for (const key in label_data) {
                    if (label_data.hasOwnProperty(key)) {
                        // Add the key (family name) to the beginning of the values
                        if (Family_exclude.includes(key)) {
                            values.push("");
                            // Add the associated values
                            values = values.concat(label_data[key]);
                        } else {
                            values.push("");
                            values.push(key);
                            Family_list.push(key)
                            // Add the associated values
                            values = values.concat(label_data[key]);
                        }
                    }
                }
            } else {
                // If Spacing is false or there is only one key
                values = Object.values(label_data).flat(); // Flatten the data without placeholders
            }
        } else {
            if (Spacing && Object.keys(label_data).length > 1) {
                // If Spacing is true and there are multiple keys
                for (const key in label_data) {
                    if (label_data.hasOwnProperty(key)) {
                        values = values.concat(label_data[key], ""); // Add an empty string as a placeholder
                    }
                }
            } else {
                // If Spacing is false or there is only one key
                values = Object.values(label_data).flat(); // Flatten the data without placeholders
            }
        }

        // Add in Class
        if (odorant) {
            Header_list = ["O2","O1"];
            if (level === 0) {
                values.unshift("");
                values.unshift("O2");
                values.push("");
            } else if (level === 1) {
                values.unshift("");
                values.unshift("O2");
                values.push("");
            } else if (level === 2) {
                values.unshift("");
                values.unshift("O2");
                values.push("");
            } else if (level === 3) {
                values.unshift("");
                values.unshift("O1");
                values.push("");
            }

        } else {
            Header_list = ["A","B1","B2","C","F","T2","Classless"];
            if (level === 0) {
                values.unshift("");
                values.unshift("A");
                values.push("");
            } else if (level === 1) {
                values.unshift("A");
                values.push("");
            } else if (level === 2) {
                if (Spacing) {
                    values.splice(26, 0, "B2");
                    values.splice(27, 0, "");
                    values.push("")

                } else {
                    values.splice(15, 0, "B2");
                }
                values.unshift("B1");
            } else if (level === 3) {
                values.unshift("C");
                values.splice(33, 0, "F");
                values.splice(33, 0, "");
                values.push("")
            } else if (level === 4) {
                values.unshift("T2");
                values.splice(28, 0, "Classless");
                values.splice(28, 0, "");
                values.splice(30, 0, "");
                values.push("")
            } else if (level === 5) {
                values.unshift("CLASS F");
            } else if (level === 6) {
                values.unshift("CLASSLESS");
            }
        }
        // Header_list = ["CLASS A","A CONT.","A ORPHANS","CLASS B1","CLASS B2","CLASS C","CLASS F","CLASS T2","CLASSLESS"];

        function calculatePositionAndAngle(index, total, values, Header_list, Family_list, isSplit) {
            // Get the text value at the current index
            const text_value = values[index];

            // Check if the text value is in the Header_list
            const isInHeaderList = Header_list.includes(text_value)  && text_value !== "Classless";

            // Check if the text value is in the Family_list
            const isInFamilyList = Family_list.includes(text_value);

            // Offset the angle calculation by -90 degrees (or -π/2 radians) to start at 12 o'clock
            const angle = -((index / total) * 2 * Math.PI) + (Math.PI / 2);

            // If isSplit is true, treat it as a Family_list item (or handle it in a special way)
            const adjustedRadius = isSplit ? (GPCRome_radius - 8) : text_value === "Classless" ? (GPCRome_radius - 10) : (isInFamilyList ? (GPCRome_radius - 20) : (isInHeaderList ? (GPCRome_radius + 18) : (GPCRome_radius + label_offset)));

            // Position on the GPCRome's border with or without label offset
            const x = width / 2 + Math.cos(angle) * adjustedRadius;
            const y = height / 2 - Math.sin(angle) * adjustedRadius;

            // If it's a header, set the rotation to 0, otherwise calculate the outward-facing rotation
            let rotation;
            if (isInHeaderList) {
                rotation = 0;
            } else if (text_value === "Classless") {
                rotation = 0;
            } else {
                rotation = -(angle * 180 / Math.PI);
            }

            return { x, y, rotation };
        }

       // Bind data and append text elements for the specific GPCRome
        svg.selectAll(`.GPCRome-text-${level}`)
        .data(values)
        .enter()
        .append("text")
        .attr("class", (d) => {
            let baseClass = `GPCRome-text-${level}`;
            // Add highlight class if the label is in the Header_list
            if (Header_list.includes(d)) {
                baseClass += ` GPCRome-text-${level}-highlight`;
            }
            // Add a family-specific class if the label is in the Family_list
            if (Family_list.includes(d)) {
                baseClass += " GPCRome-family-label";  // Add this class for family labels
            }
            return baseClass;
        })
        .attr("x", (d, i) => {
            const pos = calculatePositionAndAngle(i, values.length, values, Header_list, Family_list, false);
            return pos.x;
        })
        .attr("y", (d, i) => {
            const pos = calculatePositionAndAngle(i, values.length, values, Header_list, Family_list, false);
            return pos.y;
        })
        .attr("text-anchor", (d, i) => {
            // Center the text for headers, and handle normal text alignment for others
            if (Header_list.includes(d) && d !== "Classless" && d !== "A cont.") {
                return "middle";  // Horizontally center the headers
            }
            const angle = (i / values.length) * 360 - 90;
            return (angle >= -90 && angle < 90) ? "start" : "end";
        })
        .attr("dominant-baseline", "middle")
        .attr("dy", (d) => Header_list.includes(d) ? "0.1em" : "0.05em")  // Adjust 'dy' as needed
        .attr("transform", (d, i) => {
            const pos = calculatePositionAndAngle(i, values.length, values, Header_list, Family_list, false);

            // Calculate the angle and determine the text's side (right or left)
            const angle = (i / values.length) * 360 - 90;

            // Rotation logic
            let rotation;
            if (Header_list.includes(d)) {
                // Headers have no rotation (0 degrees)
                rotation = 0;
            } else {
                // For non-headers, flip the text on the left-hand side by 180 degrees
                rotation = angle >= -90 && angle < 90 ? 0 : 180;
            }

            // Apply the rotation and positioning
            return `rotate(${pos.rotation + rotation}, ${pos.x}, ${pos.y})`;
        })
        .html(d => GPCRome_formatTextWithHTML(d, Family_list, level))
        .style("font-size", d => Header_list.includes(d) ? Class_fontsize : Receptor_and_family_fontsize)
        .style("font-family", Font_family)
        .attr("font-weight", d => Header_list.includes(d) || Family_list.includes(d) ? "950" : "normal")
        .style("fill", d => Header_list.includes(d) ? "Black" : "black");

        // Function to process the formattedText
        function formatText(text) {
            if (text.includes("Odorant")) {
                // Split the text by the word "Odorant" and trim any leading/trailing spaces
                let result = text.split("Odorant").pop().trim();

                // Capitalize the first letter and ensure the rest is lowercase
                return result.charAt(0).toUpperCase() + result.slice(1).toLowerCase();
            }
            return text;  // Return the text unchanged if it doesn't contain "Odorant"
        }

        // After drawing all the elements, adjust the y-position for all family labels
        // Adjust the y-position for all family labels based on the midpoint between current and previous positions
        svg.selectAll(".GPCRome-family-label")
            .each(function(d) {
                const textElement = d3v4.select(this);

                // Find the index of the family label within the full values array
                const index = values.indexOf(d);  // This gets the actual index of the current family label in the `values` array

                if (index !== -1 && index > 0) {  // Ensure the index is valid and not the first item (since we need index - 1)

                    const totalItems = values.length; // Total number of items in the current GPCRome

                    // Determine if the text anchor should be "start" or "end"
                    const angle = (index / totalItems) * 360 - 90;  // Calculate the angle based on the index
                    const additionalRotation = angle >= -90 && angle < 90 ? 0 : 180;  // Conditional rotation adjustment

                    // Format the text before checking the length
                    const formattedText = GPCRome_formatTextWithHTML(d, Family_list);

                    // Remove the existing text element before appending the split elements
                    textElement.remove();

                    // fontsize
                    let family_fontsize = Receptor_and_family_fontsize;

                    // Check if the formatted text is longer than 10 characters (or any desired length)
                    if (formattedText.length > 18) {
                        let splitIndex;
                        if (formattedText.includes("-")) {
                            // If the text contains a "-", split after the "-"
                            splitIndex = formattedText.indexOf("-",3) + 1;
                        } else {
                            // Otherwise, split at the nearest space
                            splitIndex = formattedText.lastIndexOf(" ", formattedText.length-1);
                        }
                        const firstPart = formattedText.substring(0, splitIndex);  // First part
                        const secondPart = formattedText.substring(splitIndex);  // Second part

                        // Get the current and previous positions using calculatePositionAndAngle with the isSplit flag
                        const currentPos = calculatePositionAndAngle(index, totalItems, values, Header_list, Family_list, true);
                        const prevPos = calculatePositionAndAngle(index - 1, totalItems, values, Header_list, Family_list, true);

                        off_set = level+1

                        if (angle >= -90 && angle < 90) {
                            // Right-hand side: use prevPos for the first part and currentPos for the second part

                            // Append the first part of the text (using prevPos)
                            svg.append("text")
                                .attr("x", prevPos.x)
                                .attr("y", prevPos.y+off_set)
                                .attr("text-anchor", "start")
                                .attr("dominant-baseline", "middle")
                                .attr("transform", `rotate(${prevPos.rotation + additionalRotation}, ${prevPos.x}, ${prevPos.y})`)
                                .attr("class", "GPCRome-family-label-split")
                                .text(firstPart)
                                // .style("font-weight", "bold")
                                .style("font-family", Font_family)
                                .style("font-size",family_fontsize);



                            // Append the second part of the text (using currentPos)
                            svg.append("text")
                                .attr("x", currentPos.x)
                                .attr("y", currentPos.y-off_set)
                                .attr("text-anchor", "start")
                                .attr("dominant-baseline", "middle")
                                .attr("transform", `rotate(${currentPos.rotation + additionalRotation}, ${currentPos.x}, ${currentPos.y})`)
                                .attr("class", "GPCRome-family-label-split")
                                .text(secondPart)
                                // .style("font-weight", "bold")
                                .style("font-family", Font_family)
                                .style("font-size", family_fontsize);

                        } else {
                            // Left-hand side: use currentPos for the first part and prevPos for the second part

                            // Append the first part of the text (using currentPos)
                            svg.append("text")
                                .attr("x", currentPos.x)
                                .attr("y", currentPos.y+off_set)
                                .attr("text-anchor", "end")
                                .attr("dominant-baseline", "middle")
                                .attr("transform", `rotate(${currentPos.rotation + additionalRotation}, ${currentPos.x}, ${currentPos.y})`)
                                .attr("class", "GPCRome-family-label-split")
                                .text(firstPart)
                                // .style("font-weight", "bold")
                                .style("font-family", Font_family)
                                .style("font-size",family_fontsize);

                            // Append the second part of the text (using prevPos)
                            svg.append("text")
                                .attr("x", prevPos.x)
                                .attr("y", prevPos.y-off_set)
                                .attr("text-anchor", "end")
                                .attr("dominant-baseline", "middle")
                                .attr("transform", `rotate(${prevPos.rotation + additionalRotation}, ${prevPos.x}, ${prevPos.y})`)
                                .attr("class", "GPCRome-family-label-split")
                                .text(secondPart)
                                // .style("font-weight", "bold")
                                .style("font-family", Font_family)
                                .style("font-size",family_fontsize);
                        }

                    } else {
                        // If the formatted text is shorter than 10 characters, handle it normally

                        // Get the current and previous positions without splitting (isSplit = false)
                        const currentPos = calculatePositionAndAngle(index, totalItems, values, Header_list, Family_list, false);
                        const prevPos = calculatePositionAndAngle(index - 1, totalItems, values, Header_list, Family_list, false);

                        const midX = (currentPos.x + prevPos.x) / 2;
                        const midY = (currentPos.y + prevPos.y) / 2;
                        const midRotation = (currentPos.rotation + prevPos.rotation) / 2;

                        // Append the formatted text in the middle position
                        svg.append("text")
                            .attr("x", midX)
                            .attr("y", midY)
                            .attr("dominant-baseline", "middle")
                            .attr("text-anchor", (angle >= -90 && angle < 90) ? "start" : "end")
                            .attr("transform", `rotate(${midRotation + additionalRotation}, ${midX}, ${midY})`)
                            .attr("class", "GPCRome-family-label")
                            .text(formatText(formattedText))
                            // .style("font-weight", "bold")
                            .style("font-family", Font_family)
                            .style("font-size",family_fontsize);
                    }
                }
            });

        // Define color scale for continuous data
        let colorScale;

        if (GPCRomes_styling.data_color_complexity === 'One') {
        // White to Max (One color)
        colorScale = d3v4.scaleLinear()
            .domain([GPCRomes_styling.minValue, GPCRomes_styling.maxValue])  // Only two points in the domain
            .range(['#FFFFFF', GPCRomes_styling.colorEnd]);  // White to Max color

        } else if (GPCRomes_styling.data_color_complexity === 'Two') {
        // Min to Max (Two colors)
        colorScale = d3v4.scaleLinear()
            .domain([GPCRomes_styling.minValue, GPCRomes_styling.maxValue])  // Min to Max in the domain
            .range([GPCRomes_styling.colorStart, GPCRomes_styling.colorEnd]);  // Min to Max color in range

        } else if (GPCRomes_styling.data_color_complexity === 'Three') {
        // Min to White to Max (Three colors)
        colorScale = d3v4.scaleLinear()
            .domain([GPCRomes_styling.minValue, GPCRomes_styling.avg_value, GPCRomes_styling.maxValue])  // Min, Avg, Max in the domain
            .range([GPCRomes_styling.colorStart, '#FFFFFF', GPCRomes_styling.colorEnd]);  // Min to White to Max in range
}

        // Add large hollow pie chart for the entire level
        const arcGenerator = d3v4.arc()
            .innerRadius(GPCRome_radius - 7)  // Adjust to control the hollow center size
            .outerRadius(GPCRome_radius)  // Adjust to control the thickness of the pie
            // .padAngle(level === 4 ? 0.3 : 0); // Apply padding only if level is 4

        const pieGenerator = d3v4.pie()
            .sort(null)
            .value(1)  // Create equal slices for each value
            .startAngle(-Math.PI / values.length)  // Offset to move the slices left by half their size
            .endAngle(2 * Math.PI - Math.PI / values.length);  // Correct end angle for full circle

        const pieData = pieGenerator(values);

        svg.selectAll(`.large-hollow-pie-${level}`)
            .data(pieData)
            .enter()
            .append("path")
            .attr("class", `large-hollow-pie-${level}`)
            .attr("d", arcGenerator)
            .attr("transform", `translate(${width / 2}, ${height / 2})`)
            .style("fill", (d) => {
                const value = fill_data[d.data]?.Value1;
                if (datatype === "Continuous") {
                    // Use color scale for continuous data
                    const numericValue = parseFloat(value);

                    if (numericValue === 0) {
                        return "white";  // Return "white" if the value is 0
                    }

                    return !isNaN(numericValue) ? colorScale(numericValue) : "none";
                } else if (datatype === "Discrete") {
                    // Handle discrete data or default case
                    let discrete_color_min;
                    let discrete_color_max;
                    if (GPCRomes_styling.data_color_complexity === 'One') {
                        discrete_color_min = '#FFFFFF';
                        discrete_color_max = GPCRomes_styling.colorEnd;
                    } else {
                        discrete_color_min = GPCRomes_styling.colorStart;
                        discrete_color_max = GPCRomes_styling.colorEnd;
                    }
                    if (value === "Yes") return discrete_color_max;
                    if (value === "No") return discrete_color_min;
                    return "none";  // Make the slice invisible if the value is ""
                // Expand this section for handling specific coverage pages with colors
                } else if (datatype === "Structure") {
                    // Handle discrete data or default case
                    if (value === "Active") return "green";
                    if (value === "Inactive") return "red";
                    if (value === "Both") return "blue";
                    if (value === "empty") return "white";
                    return "none";  // Make the slice invisible if the value is ""
                } else if (datatype === "Arrestin") {
                    // Handle discrete data or default case
                    if (value === "ARRB1") return "orange";
                    if (value === "ARRB2") return "purple";
                    if (value === "ARRC") return "aquamarine";
                    if (value === "ARRS") return "cornflowerblue";
                    if (value === "empty") return "white";
                    return "none";  // Make the slice invisible if the value is ""
                } else if (datatype === "NRDD") {
                    // Handle discrete data or default case
                    if (value === 1) return "#F5BCBF";
                    if (value === 2) return "#F17270";
                    if (value === 3) return "#DD2628";
                    if (value === 4) return "#2C87C8";
                    if (value === 5) return "#D3D3D3";
                    if (value === 6) return "#A3D9C8";
                    if (value === 7) return "White";
                    return "none";  // Make the slice invisible if the value is ""
                }
            })
            .style("stroke", (d) => {
                const value = fill_data[d.data]?.Value1;
                if (value === "Yes" || value === "No" || !isNaN(value)) return "black";
                if (value === "Active" || value === "Inactive" || value === "Both" || value === "empty") return "black";
                if (value === "ARRB1" || value === "ARRB2" || value === "ARRC" || value === "ARRS" || value === "empty") return "black";
                return "none";  // Remove the stroke if the value is ""
            })
            .style("stroke-width", (d) => {
                const value = fill_data[d.data]?.Value1;
                return value === "" ? 0 : 0.5;  // Set stroke-width to 0 if the value is an empty string
            });
    }
    // Add padding and scale
    const padding = 10;  // Adjust padding value as needed (50px for this example)
    const originalWidth = +svg.attr("width");
    const originalHeight = +svg.attr("height");

    // Adjust the viewBox to add padding
    svg.attr("viewBox", `-${padding} -${padding} ${originalWidth + 2 * padding} ${originalHeight + 2 * padding}`);

    // Apply scaling to the content (e.g., 95% of the original size)
    svg.attr("transform", "scale(0.95)")
    .attr("transform-origin", "center");
}
