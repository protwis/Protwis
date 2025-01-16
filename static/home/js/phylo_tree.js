function draw_tree(data, options) {

    var branches = {};
    var branch_offset = 0;
    var thickness = options.depth + 1;
    for (var key in options.branch_length) {
        if (key == options.depth) { continue };
        if (options.label_free.includes(parseInt(key))) {
            branch_offset = branch_offset + 10;
        } else {
            if (options.branch_trunc != 0) {
                branch_offset = branch_offset + 2*options.branch_trunc + 10;
            } else {
                branch_offset = branch_offset + string_pixlen(options.branch_length[key], key);
            }
        }
        branches[key] = branch_offset;
    }
    branches[options.depth] = branch_offset + options.leaf_offset;

    var diameter = 2 * branches[options.depth] + 100;

    var tree = d3.layout.tree()
        .size([360, diameter / 2])
        .separation(function (a, b) { return (a.parent == b.parent ? 1 : 2) / a.depth; });

    var diagonal = d3.svg.diagonal.radial()
        .projection(function (d) { return [d.y, d.x / 180 * Math.PI]; });

    var svg = d3.select('#'+options.anchor).append("svg")
        .attr("width", diameter)
        .attr("height", diameter)
        .attr("id", options.anchor+"_svg")
        .attr("xmlns", "http://www.w3.org/2000/svg");

    var svg_g = svg.append("g")
        .attr("transform", "translate(" + diameter / 2 + "," + diameter / 2 + ")");

    var nodes = tree.nodes(data);

    nodes.forEach(function (d) {
        if (d.depth == 0) {
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
            if ((d.target.interactions > 0 && d.target.mutations_an > 0) || 1 == 1) { return 0.8 } //|| 1==1
            else if (d.target.interactions > 0) { return 0.5 }
            else if (d.target.mutations_an > 0) { return 0.5 }
            else { return 0.1 };
        });

    var node = svg_g.selectAll(".node")
        .data(nodes)
        .enter().append("g")
        .attr("class", "node")
        .attr("transform", function (d) { if (d.name == '') { return "rotate(" + (d.x) + ")translate(" + d.y + ")"; } else { return "rotate(" + (d.x - 90) + ")translate(" + d.y + ")"; } })
//TODO: add a check to remove circles when nothing is passed (?)
    node.filter(function (d) { return (d.depth == options.depth) })
        .filter(function (d) { return (d.value !== 3000) })
        .append("circle")
        .attr("r", function (d) { if (d.name == '') { return "0" } else { return "4.0" } })
        .style("stroke", "black")
        .style("stroke-width", ".3px")
        .style("fill", function (d) {
            if (d.color && d.depth < options.depth) { return d.color }
            else if ( d.value === 1) {
                return "FireBrick";
            }
            else if ( d.value === 10) {
                return "LightGray";
            }
            else if ( d.value === 20) {
                return "DarkGray";
            }
            else if ( d.value === 30) {
                return "Gray";
            }
            else if ( d.value === 40) {
                return "Black";
            }
            else if (d.value === 100) {
                return 'LightGray';
            }
            else if (d.value === 500) {
                return 'DarkGray';
            }
            else if (d.value === 1000) {
                return 'Gray';
            }
            else if (d.value === 2000) {
                return 'Black';
            }
            else { return "White" };
        })
        .style("opacity", .99);

    node.filter(function (d) { return (d.depth == options.depth) })
        .attr("id", function (d) { if (d.name == '') { return "innerNode" } else { return 'X'+d.name.toUpperCase() } });

    node.append("text")
        .attr("dy", ".31em")
        .attr("name", function (d) { if (d.name == '') { return "branch" } else { return d.name } })
        .attr("text-anchor", function (d) {
            if (d.depth == options.depth ) {
              return d.x < 180 ? "start" : "end";
            } else {
              return d.x < 180 ? "end" : "start";
            }
        })
        .attr("transform", function (d) {
            if (d.depth == options.depth) {
                return d.x < 180 ? "translate(7)" : "rotate(180)translate(-7)";
            } else {
                return d.x < 180 ? "translate(-12)" : "rotate(180)translate(12)";
            }
        })
        .text(function (d) {
            if (d.depth == options.depth) {
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

        .style("font-size", function (d) { if (d.depth < 2) { return "14px" } else if (d.depth == 2) { return "12px" } else { return "10px" } })
        .style("font-family", "Palatino")
        .style("fill", function (d) {
            if (d.color) { return "#111" }
            else { return "#222" };
        }).call(getBB);
    node.filter(function (d) { return (d.depth != options.depth) }).insert("rect", "text")
        .attr("x", function (d) { return d.x < 180 ? d.bbox.x - 12 : d.bbox.x - d.bbox.width - 12; })
        .attr("y", function (d) { return d.bbox.y })
        .attr("width", function (d) { return d.bbox.width })
        .attr("height", function (d) { return d.bbox.height })
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
        if (depth < 2) {
            ctx.font = "20px Palatino"
        } else if (depth == 2) {
            ctx.font = "14px Palatino"
        } else {
            ctx.font = "12px Palatino"
        }
        return parseInt(ctx.measureText(text).width) + 40;
    }

    function getBB(selection) {
        selection.each(function (d) { d.bbox = this.getBBox(); })
    }

    function wrap(text, width) {
        if (width == 0) {
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
            while (word = words.pop()) {
                line.push(word);
                tspan.text(line.join(" "));
                if (tspan.node().getComputedTextLength() > width) {
                    line.pop();
                    tspan.text(line.join(" "));
                    line = [word];
                    tspan = text.append("tspan").attr("x", 0).attr("y", y).attr("dy", ++lineNumber * lineHeight + dy + "em").text(word);
                }
            }
        });
    }
}


/**
* changeLeavesLabels
*
* Function designed to change label names of phylo tree
*
* @location {string} location - svg in which to draw outer circles
* @value {string} value - either IUPHAR or UniProt for label names
* @dict {dictionary} dict - the translation dictionary for label names (from IUPHAR to UniProt and vice versa)
*/

function changeLeavesLabels(location, value, dict){
  // Initialize leaf node length
  maxLeafNodeLenght = 0;
  // Find longest label
  gNodes = d3.select('#'+location).selectAll('g');
  gNodes.each(function(d) {
    if (d3.select(this).attr("id") !== null) {
      name = d3.select(this).attr("id").substring(1);
      labelName = dict[name][0];
      // replaces labels derived from view
      labelName = labelName.replace("-adrenoceptor", '');
      labelName = labelName.replace(" receptor-", '-');
      labelName = labelName.replace("<sub>", '</tspan><tspan baseline-shift = "sub">');
      labelName = labelName.replace("</sub>", '</tspan><tspan>');
      labelName = labelName.replace("<i>", '</tspan><tspan font-style = "italic">');
      labelName = labelName.replace("</i>", '</tspan><tspan>');
      labelName = labelName.replace("Long-wave-sensitive",'LWS');
      labelName = labelName.replace("Medium-wave-sensitive",'MWS');
      labelName = labelName.replace("Short-wave-sensitive",'SWS');
      labelName = labelName.replace("Olfactory", 'OLF');
      labelName = labelName.replace("Calcitonin -like", 'CLR');
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

/**
* DrawCircles
*
* Function designed to append data circles on the external part of phylo tree
*
* @location {string} location - svg in which to draw outer circles
* @data {Object} data - data provided by the view (json dict usually)
* @starter {integer} starter - the max length of the leaves, to start drawing the circles (calculated by changeLeavesLabels)
* @dict {dictionary} dict - the translation dictionary for color codes
* @fancy {boolean} fancy - set the option for fancy circles
* @clean {boolean} clean - remove the inner circles in the plot
*/

function DrawCircles(location, data, starter, dict, fancy=false, clean=true, fixed_total=false, gradient=true){

    // const pSBC=(p,c0,c1,l)=>{
    //     let r,g,b,P,f,t,h,i=parseInt,m=Math.round,a=typeof(c1)=="string";
    //     if(typeof(p)!="number"||p<-1||p>1||typeof(c0)!="string"||(c0[0]!='r'&&c0[0]!='#')||(c1&&!a))return null;
    //     if(!this.pSBCr)this.pSBCr=(d)=>{
    //         let n=d.length,x={};
    //         if(n>9){
    //             [r,g,b,a]=d=d.split(","),n=d.length;
    //             if(n<3||n>4)return null;
    //             x.r=i(r[3]=="a"?r.slice(5):r.slice(4)),x.g=i(g),x.b=i(b),x.a=a?parseFloat(a):-1
    //         }else{
    //             if(n==8||n==6||n<4)return null;
    //             if(n<6)d="#"+d[1]+d[1]+d[2]+d[2]+d[3]+d[3]+(n>4?d[4]+d[4]:"");
    //             d=i(d.slice(1),16);
    //             if(n==9||n==5)x.r=d>>24&255,x.g=d>>16&255,x.b=d>>8&255,x.a=m((d&255)/0.255)/1000;
    //             else x.r=d>>16,x.g=d>>8&255,x.b=d&255,x.a=-1
    //         }return x};
    //     h=c0.length>9,h=a?c1.length>9?true:c1=="c"?!h:false:h,f=this.pSBCr(c0),P=p<0,t=c1&&c1!="c"?this.pSBCr(c1):P?{r:0,g:0,b:0,a:-1}:{r:255,g:255,b:255,a:-1},p=P?p*-1:p,P=1-p;
    //     if(!f||!t)return null;
    //     if(l)r=m(P*f.r+p*t.r),g=m(P*f.g+p*t.g),b=m(P*f.b+p*t.b);
    //     else r=m((P*f.r**2+p*t.r**2)**0.5),g=m((P*f.g**2+p*t.g**2)**0.5),b=m((P*f.b**2+p*t.b**2)**0.5);
    //     a=f.a,t=t.a,f=a>=0||t>=0,a=f?a<0?t:t<0?a:a*P+t*p:0;
    //     if(h)return"rgb"+(f?"a(":"(")+r+","+g+","+b+(f?","+m(a*1000)/1000:"")+")";
    //     else return"#"+(4294967296+r*16777216+g*65536+b*256+(f?m(a*255):0)).toString(16).slice(1,f?undefined:-2)
    // }

    // Function to calculate the color shade based on percentage
    function calculateColorShade(percentage, colorCode) {
      // Parse the color code to RGB
      const r = parseInt(colorCode.slice(1, 3), 16);
      const g = parseInt(colorCode.slice(3, 5), 16);
      const b = parseInt(colorCode.slice(5, 7), 16);

      // Calculate the new RGB values based on the percentage
      const newR = Math.round(r + (255 - r) * (1-percentage));
      const newG = Math.round(g + (255 - g) * (1-percentage));
      const newB = Math.round(b + (255 - b) * (1-percentage));

      // Convert the RGB values back to a hexadecimal color with two digits per channel
      const newColorCode = `#${newR.toString(16).padStart(2, '0')}${newG.toString(16).padStart(2, '0')}${newB.toString(16).padStart(2, '0')}`;

      return newColorCode;
    }

    function componentToHex(c)
    {
        var hex = c.toString(16);
        return hex.length == 1 ? "0" + hex : hex;
    }

    var svg = d3.select('#'+location);
    var node = svg.selectAll(".node");
    if (clean === true) {
      node.selectAll("circle").remove();
    }
    if (fancy === false) {
      var spacer = 8;
        for (var x in data){
          for (var unit in dict){
              // if (data[x].constructor == Object) {
              //     var index = Object.keys(data[x]).indexOf(unit);
              // }
              // else {
              var index = data[x].indexOf(unit);
              // }
              if (index >= 0) {
                  // variable to set the location of the different circle drawing
                  multiply = 1+Object.keys(dict).indexOf(unit);
                  var leafwithname = svg.selectAll('g[id=X'+x+']')
                      .append("circle")
                      .attr("r", 3.25)
                      .style("fill", dict[unit])
                      .attr("transform", "translate(" + (Math.ceil(starter) + multiply*spacer) + ",0)");
              }
          }
        }
    } else {
      var spacer = 10;
      node.selectAll("ellipse").remove();
      for (var x in data){
        keys = Object.keys(data[x]);
        tot = Object.keys(data[x]).reduce((a, b) => a + b, 0);
        var sum = 0;
        for( var el in data[x] ) {
          if( data[x].hasOwnProperty( el ) ) {
            sum += parseFloat( data[x][el] );
          }
        }
        for (var unit in dict){
          if (keys.indexOf(unit)>= 0) {
            if (gradient) {
                // Caculate percentage with fixed total
                if (fixed_total) {
                    total = fixed_total;
                }
                else {
                    total = Object.values(data[x]).reduce((acc, val) => acc + val, 0);
                }
                value = data[x][unit];
                percentage = (value / total) * 100;
                // Overwrite percentage to have empty circles for entries with data length = 1
                if (fixed_total && Object.keys(data[x]).length===1) {
                    percentage = 0;
                }
                // color= pSBC(percentage, dict[unit]);
                color = calculateColorShade(percentage / 100, dict[unit]);
            }
            else {
                color = dict[unit];
            }
            multiply = 1+Object.keys(dict).indexOf(unit);
            leaf = svg.selectAll('g[id=X'+x+']');
            var leafwithname = svg.selectAll('g[id=X'+x+']')
                .append("ellipse")
                .attr("rx", 4.25)
                .attr("ry", 4.25)
                .style("stroke", dict[unit])
                .style("stroke-width", 1.2)
                .style("fill", color)
                .attr("transform", "translate(" + (Math.ceil(starter) + multiply*spacer) + ",0)");
            }
        }
      }
    }
  }

/**
* FancyCircles
*
* Function designed to append fancy data circles on the external part of phylo tree with
*
* @location {string} location - svg in which to draw outer circles
* @data {Object} data - data provided by the view (json dict usually)
* @starter {integer} starter - the max length of the leaves, to start drawing the circles (calculated by changeLeavesLabels)
* @dict {dictionary} dict - the translation dictionary for color codes
*/

function FancyCircles(location, data, starter, dict){
    var spacer = 8;
    var svg = d3.select('#'+location);
    var node = svg.selectAll(".node");
    for (var x in data){
      keys = Object.keys(data[x]);
      max = Object.keys(data[x]).reduce(function(a, b){ return data[x][a] > data[x][b] ? a : b });
      for (var unit in dict){
        if (keys.indexOf(unit)>= 0) {
          // variable to set the location of the different circle drawing
          multiply = 1+Object.keys(dict).indexOf(unit);
          if (unit === max) {
            var leafwithname = svg.selectAll('g[id=X'+x+']')
                .append("circle")
                .attr("r", 3.25)
                .style("stroke", dict[unit])
                .style("stroke-width", 1)
                .style("fill", dict[unit])
                .attr("transform", "translate(" + (Math.ceil(starter) + multiply*spacer) + ",0)");
          }
          else {
            var leafwithname = svg.selectAll('g[id=X'+x+']')
                .append("circle")
                .attr("r", 3.25)
                .style("stroke", dict[unit])
                .style("stroke-width", 1)
                .style("fill", "white")
                .attr("transform", "translate(" + (Math.ceil(starter) + multiply*spacer) + ",0)");
              }
          }
      }
    }
  }

/**
* draw_cluster
*
* Function designed to draw horizontal cluster
*
*
* @data {Object} data - data provided by the view (json dict usually)
* @options {dictionary} options - options provided for the sake of depth, data info and text anchor
* @trim {boolean} trim - boolean value to select trimming empty leaves
*/

function draw_cluster(data, options, trim=true) {

    if (trim==true) {
      for (var family in data['children']) {
        for (var child in data['children'][family]['children']){
          var i = data['children'][family]['children'][child]['children'].length;
          while (i-- ){
            if (data['children'][family]['children'][child]['children'][i]['value'] == 0) {
              data['children'][family]['children'][child]['children'].splice(i, 1);
            }
          }
        }
        var k = data['children'][family]['children'].length;
        while (k-- ){
          if (data['children'][family]['children'][k]['children'].length == 0) {
            data['children'][family]['children'].splice(k, 1);
          }
        }
      }
      var j = data['children'].length;
      while (j-- ){
        if (data['children'][j]['children'].length == 0) {
          data['children'].splice(j, 1);
        }
      }
    }

    var branches = {};
    var branch_offset = 0;
    for (var key in options.branch_length) {
        if (key == options.depth) { continue };
        if (options.label_free.includes(parseInt(key))) {
            branch_offset = branch_offset + 10;
        } else {
            if (options.branch_trunc != 0) {
                branch_offset = branch_offset + 2*options.branch_trunc + 10;
            } else {
                branch_offset = branch_offset + string_pixlen(options.branch_length[key], key);
            }
        }
        branches[key] = branch_offset;
    }
    branches[options.depth] = branch_offset + options.leaf_offset;

    var dimension = 2 * branches[options.depth] + 100;

    var cluster = d3.layout.cluster()
        .size([dimension/3, (dimension/3) - 100])
        .separation(function (a, b) { return (a.parent == b.parent ? 1 : 2) / a.depth; });

    var diagonal = d3.svg.diagonal()
        .projection(function (d) {
          return [d.y, d.x];
        });

    var svg = d3.select('#'+options.anchor).append("svg")
        .attr("width", dimension)
        .attr("height", dimension/2)
        .attr("id", options.anchor+"_svg")
        .attr("transform", "translate(0," + (dimension/12) + ")");



    var nodes = cluster.nodes(data);

    nodes.forEach(function (d) {
        if (d.depth == 0) {
            d.y = 0
        } else {
            d.y = branches[d.depth]
        }
    });

    var links = cluster.links(nodes);

    var link = svg.selectAll(".link")
        .data(links)
        .enter()
        .append("path")
        .each(function (d) { d.target.linkNode = this; })
        .attr("class", "link_continuous")
        .attr("d", diagonal)
        .style("stroke", function (d) { return d.target.color; })
        .style("stroke-width", function (d) { if (d.target.depth > 0) { return 4 - d.target.depth; } else { return 0; } })
        .style("fill-opacity", 0);

    var node = svg.selectAll(".node")
        .data(nodes)
        .enter()
        .append("g")
        .attr("class", "node")
        .attr("transform", function (d) {
      return "translate(" + d.y + "," + d.x + ")";
        })

    node.filter(function (d) { return (d.depth == options.depth) })
        .append("circle")
        .attr("r", "4.0")
        .style("stroke", "black")
        .style("stroke-width", ".3px")
        .style("fill", function (d) {
            if (d.color && d.depth < options.depth) { return d.color }
            else if ( d.value === 10) {
                return "LightGray";
            }
            else if ( d.value === 20) {
                return "DarkGray";
            }
            else if ( d.value === 30) {
                return "Gray";
            }
            else if ( d.value === 40) {
                return "Black";
            }
            else { return "White" };
        })
        .style("opacity", .99);

    node.filter(function (d) { return (d.depth == options.depth) })
        .attr("id", function (d) { if (d.name == '') { return "innerNode" } else { return 'X'+d.name.toUpperCase() } });

    node.append("text")
        .attr("dy", "3")
        .attr("dx", function(d) { return d.children ? -8 : 8; })
        .attr("name", function (d) { if (d.name == '') { return "branch" } else { return d.name } })
        .attr("text-anchor", function (d) {
              return d.children ? "end" : "start";
        })
        .text(function (d) {
            if (d.depth == options.depth) {
                return d.name.toUpperCase();
            } else if (options.label_free.includes(d.depth)) {
                return "";
            } else if (d.depth > 0) {
                return d.name;
            } else {
                return "";
            }
        })
        // .call(wrap, options.branch_trunc)
        .style("font-size", "12px")
        .style("font-family", "Palatino")
        .style("fill", function (d) {
            if (d.color) { return "#111" }
            else { return "#222" };
        }).call(getBB);
    node.filter(function (d) { return (d.depth != options.depth) }).insert("rect", "text")
        .attr("x", function (d) { return d.bbox.x })
        .attr("y", function (d) { return d.bbox.y })
        .attr("width", function (d) { return d.bbox.width })
        .attr("height", function (d) { return d.bbox.height })
        .style("fill", "#FFF");

    function string_pixlen(text, depth) {
        var canvas = document.createElement('canvas');
        var ctx = canvas.getContext("2d");
        if (depth < 2) {
            ctx.font = "22px Palatino"
        } else if (depth == 2) {
            ctx.font = "16px Palatino"
        } else {
            ctx.font = "14px Palatino"
        }
        return parseInt(ctx.measureText(text).width) - 100;
    }

    function getBB(selection) {
        selection.each(function (d) { d.bbox = this.getBBox(); })
    }

    function wrap(text, width) {
        if (width == 0) {
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
                dy = parseFloat(text.attr("dy"));
                // tspan = text.text(null).append("tspan").attr("x", 0).attr("y", y).attr("dy", dy + "em");
            while (word = words.pop()) {
                line.push(word);
                tspan.text(line.join(" "));
                if (tspan.node().getComputedTextLength() > width) {
                    line.pop();
                    tspan.text(line.join(" "));
                    line = [word];
                    tspan = text.append("tspan").attr("x", 0).attr("y", y).attr("dy", ++lineNumber * lineHeight + dy + "em").text(word);
                }
            }
        });
    }

}



// Heatmap (circles_bal, heatmap_data, master_dict, buttonObjects);
function draw_heatmap(square_data, data, bible, options, location, element_id, legend_label, trim=true) {
  // set the dimensions and margins of the graph

  var support = (function () {
  	if (!window.DOMParser) return false;
  	var parser = new DOMParser();
  	try {
  		parser.parseFromString('x', 'text/html');
  	} catch(err) {
  		return false;
  	}
  	return true;
  })();

  function oddOrEven(x) {
    return ( x & 1 ) ? "odd" : "even";
  }

  var textToHTML = function (str) {

  	// check for DOMParser support
  	if (support) {
  		var parser = new DOMParser();
  		var doc = parser.parseFromString(str, 'text/html');
  		return doc.body.innerHTML;
  	}

  	// Otherwise, create div and append HTML
  	var dom = document.createElement('div');
  	dom.innerHTML = str;
  	return dom;

  };

  var counting = function (var1, var2, var3) {
    for (var i=0; i < var1.length; i++){
      id = var1[i].replace(' (','').replace(')','');
      var2[id] = (i + 1)*var3;
    }
    return var2;
  };

  var margin = {top: 30, right: 30, bottom: 30, left: 30};

  // Labels of row and columns
  var rows = [];
  var html_converter = [];
  var cols = [];
  var colorCols = ['color'];
  var chartData = [];
  var colorData = [];
  var clone = [];
  if (typeof structuredClone === "function") {
    clone = structuredClone(data);
  } else {
    clone = JSON.parse(JSON.stringify(data));
  }
  var gpcrClasses = [];
  var gpcrLigandType = [];
  var gpcrReceptorFamily = [];
  var gpcrName = [];
  var sorterDict = {};
  var filtered_square_data = {};
  var highest_value = 0;
  var rows_to_uniprot = {};

  if (trim==true) {
    for (var prot_class in clone){
      for (var entry in clone[prot_class]) {
        for (var family in clone[prot_class][entry]['children']) {
          var i = clone[prot_class][entry]['children'][family]['children'].length;
          while (i-- ){
            if (clone[prot_class][entry]['children'][family]['children'][i]['value'] == 0) {
              clone[prot_class][entry]['children'][family]['children'].splice(i, 1);
            }
          }
        }
        var k = clone[prot_class][entry]['children'].length;
        while (k-- ){
          if (clone[prot_class][entry]['children'][k]['children'].length == 0) {
            clone[prot_class][entry]['children'].splice(k, 1);
          }
        }
      }
      var l = clone[prot_class].length;
      while (l-- ){
        if (clone[prot_class][l]['children'].length == 0) {
          clone[prot_class].splice(l, 1);
        }
      }
    }
    for (var prot_class in clone) {
      gpcrClasses.push(prot_class);
      if (options['classClick'].includes(prot_class)){
        for (var entry in clone[prot_class]){
          gpcrLigandType.push(clone[prot_class][entry]['name']);
          for (var family in clone[prot_class][entry]['children']) {
            gpcrReceptorFamily.push(clone[prot_class][entry]['children'][family]['name']);
            for (var i in clone[prot_class][entry]['children'][family]['children']){
              tmp = clone[prot_class][entry]['children'][family]['children'][i]['name'].toUpperCase();
              row_val = bible[clone[prot_class][entry]['children'][family]['children'][i]['name'].toUpperCase()][options['namesClick']];
              gpcrName.push(tmp);
              filtered_square_data[tmp] = square_data[tmp];
              rows_to_uniprot[row_val] = tmp;
              rows.push(row_val);
              html_converter.push(
                {
                  row: bible[clone[prot_class][entry]['children'][family]['children'][i]['name'].toUpperCase()][options['namesClick']],
                  html: textToHTML(bible[clone[prot_class][entry]['children'][family]['children'][i]['name'].toUpperCase()][options['namesClick']])
                }
              );
              // console.log(textToHTML(bible[clone[prot_class][entry]['children'][family]['children'][i]['name'].toUpperCase()][options['namesClick']]));
              colorData.push(
                { row: bible[clone[prot_class][entry]['children'][family]['children'][i]['name'].toUpperCase()][options['namesClick']],
                  col: 'color',
                  value: bible[clone[prot_class][entry]['children'][family]['children'][i]['name'].toUpperCase()][options['colorClick']]
                }
              );
            }
          }
        }
      }
    }
  } else {
    for (var prot_class in data) {
      gpcrClasses.push(prot_class);
      if (options['classClick'].includes(prot_class)){
        for (var entry in data[prot_class]){
          gpcrLigandType.push(data[prot_class][entry]['name']);
          for (var family in data[prot_class][entry]['children']) {
            gpcrReceptorFamily.push(data[prot_class][entry]['children'][family]['name']);
            for (var i in data[prot_class][entry]['children'][family]['children']){
              gpcrName.push(data[prot_class][entry]['children'][family]['children'][i]['name'].toUpperCase());
              rows.push(
                textToHTML(bible[data[prot_class][entry]['children'][family]['children'][i]['name'].toUpperCase()][options['namesClick']])
              );
              colorData.push(
                { row: bible[data[prot_class][entry]['children'][family]['children'][i]['name'].toUpperCase()][options['namesClick']],
                  col: 'color',
                  value: bible[data[prot_class][entry]['children'][family]['children'][i]['name'].toUpperCase()][options['colorClick']]
                }
              );
            }
          }
        }
      }
    }
  }

  sorterDict = counting(gpcrClasses.sort().reverse(), sorterDict, 1000000);
  sorterDict = counting(gpcrLigandType.sort().reverse(), sorterDict, 10000);
  sorterDict = counting(gpcrReceptorFamily.sort().reverse(), sorterDict, 100);
  sorterDict = counting(gpcrName.sort().reverse(), sorterDict, 1);

  var sorted_rows = [];
  var tmpdict = {};
  for (var i=0; i < rows.length; i++){
    k = 0;
    for (var id in options['sortClick']){
      name = options['sortClick'][id];
      fam = bible[rows_to_uniprot[rows[i]]][name].replace(' receptors','').replace('(','').replace(')','');
      k = k + sorterDict[fam];
    }
    tmpdict[rows[i]] = k;

  }
  // from dictionary to array
  var tmpItems = Object.keys(tmpdict).map(function(key) {
    return [key, tmpdict[key]];
  });
  // sorting the generated array based on the calculated values
  tmpItems.sort(function(first, second) {
    return second[1] - first[1];
  });
  // pushing the sorted names to the rows
  for (var i=0; i < tmpItems.length; i++){
    sorted_rows.push(tmpItems[i][0]);
  }
  // reversing the array for plotting purposes
  sorted_rows.reverse();
  // pushing data to the chart to be plotted
  for (var receptor in filtered_square_data){
    for (var subtype in filtered_square_data[receptor]){
      if (subtype != 'null'){
        cols.push(subtype);
        chartData.push({ row: bible[receptor][options['namesClick']], col: subtype, value: filtered_square_data[receptor][subtype]});
        if (filtered_square_data[receptor][subtype] > highest_value){
          highest_value = filtered_square_data[receptor][subtype];
        }
      }
    }
  }

  var cols = cols.filter(function(item, pos){
                          return cols.indexOf(item)== pos;
                          });
  cols.sort();

  var width = (35 * cols.length);
  var height = (20 * rows.length);

  // append the svg object to the body of the page
  var svg_home = d3v4.select("#" + location)
              .append("svg")
              .attr("width", width + margin.left + margin.right)
              .attr("height", height + (margin.top*5))
              .attr("transform", "translate(0,-" + margin.bottom + ")")
              .attr("id", element_id)

  var greyscale = [ 'rgb(255,255,255)', 'rgb(0,0,0)' ];

  var legend = svg_home.append('defs')
                       .append('linearGradient')
                       .attr('id', 'grad_' + element_id )
                       .attr('x1', '0%')
                       .attr('x2', '100%')
                       .attr('y1', '0%')
                       .attr('y2', '0%');

   legend.selectAll('stop')
         .data(greyscale)
         .enter()
         .append('stop')
         .style('stop-color', function(d){ return d; })
         .attr('offset', function(d,i){
                return 100 * (i / (greyscale.length - 1)) + '%';
              })


  var legend_svg = svg_home.append("g")
                    .attr("transform", "translate(0," + (height+100)  + ")");

  var color_svg = svg_home.append("g")
                    .attr("transform", "translate(" + (margin.left*1.5) + ",0)");

  var svg = svg_home.append("g")
                    .attr("transform", "translate(" + (margin.left*2.5) + ",0)");

  legend_svg.append("text")
            // .attr('x', 75)
            .attr('y', 10)
            .style("font", "14px sans-serif")
            .style("font-weight", "bold")
            .text(legend_label);

  legend_svg.append("text")
            .attr('x', 30)
            .attr('y', 30)
            .style("font", "14px sans-serif")
            .style("font-weight", "bold")
            .text('0');

  legend_svg.append('rect')
            .attr('x', 40)
            .attr('y', 15)
            .attr('width', (width*0.9) + margin.left + margin.right)
            .attr('height', 20)
            .style('fill', 'url(#grad_' + element_id + ')');

  legend_svg.append("text")
            .attr('x', width + margin.left + (margin.right*1.75))
            .attr('y', 30)
            .style("font", "14px sans-serif")
            .style("font-weight", "bold")
            .text(highest_value);

  legendLabel = legend_svg.select('text')['_groups'][0][0].getBBox().width*1.05 + 0.5 * 10;

  legend_svg.select("text")
            .attr("x", (width + margin.left + margin.right - legendLabel)/2 + margin.left)

  // Build X scales and axis:
  var x = d3v4.scaleBand()
    .range([ 0, width ])
    .domain(cols)
    .padding(0.01);

  var xColor = d3v4.scaleBand()
    .range([ 0, width ])
    .domain(colorCols)
    .padding(0.01);

  svg.append("g")
    .attr("transform", "translate(0," + height + ")")
    .attr('id', 'Xaxis')
    .call(d3v4.axisBottom(x).tickSizeOuter(0));

  // Build X scales and axis:
  var y = d3v4.scaleBand()
    .range([ height, 0 ])
    .domain(sorted_rows)
    .padding(0.01);

  color_svg.append("g")
    .attr('id', 'Yaxis')
    .call(d3v4.axisLeft(y).tickSize(0));

  // Build color scale
  var myColor = d3v4.scaleLinear()
    .range(["white", "black"])
    .domain([1,highest_value])

  //Read the data
    color_svg.selectAll()
        .data(colorData, function(d) {return d.row+':'+d.col;})
        .enter()
        .append("rect")
        .attr("x", function(d) { return x(d.col) })
        .attr("y", function(d) { return y(d.row) })
        .attr("width", "16px" )
        .attr("height", y.bandwidth() )
        .style("fill", function(d) { return d.value} )

  //Read the data
    svg.selectAll()
        .data(chartData, function(d) {return d.row+':'+d.col;})
        .enter()
        .append("rect")
        .attr("x", function(d) { return x(d.col) })
        .attr("y", function(d) { return y(d.row) })
        .attr("width", x.bandwidth() )
        .attr("height", y.bandwidth() )
        .style("fill", function(d) { return myColor(d.value)} )

    svg.select('#Xaxis')
       .attr('text-anchor', 'start')

    d3v4.selectAll("#Yaxis>.tick>text")
        .each(function(d, i){
          d3.select(this).style("font-size","1.3em");
        });

    d3v4.selectAll("#Xaxis>.tick>text")
        .each(function(d, i){
          d3.select(this).style("font-size","1.3em");
        });

    var count = 1;
    ticks = svg.select('#Xaxis').selectAll('.tick');
    ticks.each(function(d) {
      value = oddOrEven(count);
      text = d3.select(this).select('text');
      textSize = Math.floor(text[0][0].getBBox().width*1.05 + 0.5 * 10);
      text.attr("transform", null);
      text.attr("x", "-" + textSize/2 + "px");
      if(value === "even"){
        // text.attr("x", "-15");
        text.attr("y", "30");
      }
      count = count + 1;
    });

    yTicks = color_svg.select('#Yaxis').selectAll('.tick');
    yTicks.each(function(d) {
      yText = d3.select(this).select('text')[0];
      for (var i=0; i < html_converter.length; i++){
        if (html_converter[i]['row'] === yText[0].innerHTML.replaceAll("&lt;", "<").replaceAll("&gt;", ">").replaceAll("amp;","")){
          yText[0].innerHTML = html_converter[i]['html'];
          labelSize = yText[0].getBBox().width*1.05 + 0.5 * 10
          svg_home.attr("width", width + margin.left + margin.right + labelSize);
          svg.attr("transform", "translate(" + (labelSize+40) + "," + (margin.top*1.5) + ")");
          color_svg.attr("transform", "translate(" + (labelSize+20) + "," + (margin.top*1.5) + ")");
        }
      }
    });

}


function draw_model_scores(location, element_id, score, startValue, endValue, definedValue, gradientChange, inverted=false) {

    var margin = {top: 30, right: 30, bottom: 30, left: 30};
    var width = 200;
    var height = 50;

    // append the svg object to the body of the page
    var svg_home = d3v4.select("#" + location)
                .append("svg")
                .attr("width", width + margin.left + margin.left + margin.right)
                .attr("height", (margin.top))
                .attr("transform", "translate(-30,0)")
                .attr("id", element_id);

    if (inverted == true){
        var greyscale = ['rgb(0,0,0)', 'rgb(255, 255, 255)' ];
    } else {
        var greyscale = [ 'rgb(255, 255, 255)', 'rgb(0,0,0)' ];
    };

    var first_legend = svg_home.append('defs')
                         .append('linearGradient')
                         .attr('id', 'grad_' + score)
                         .attr('x1', '0%')
                         .attr('x2', '100%')
                         .attr('y1', '0%')
                         .attr('y2', '0%');

    var stopOffset = (gradientChange - startValue) / (endValue - startValue) * 100;

    first_legend.append('stop')
           .style('stop-color', greyscale[0])
           .attr('offset', '0%');
    // first_legend.append('stop')
    //        .style('stop-color', greyscale[1])
    //        .attr('offset', stopOffset + '%');
    // first_legend.append('stop')
    //        .style('stop-color', greyscale[1])
    //        .attr('offset', stopOffset + '%');
    first_legend.append('stop')
           .style('stop-color', greyscale[1])
           .attr('offset', '100%');


    var first_legend_svg = svg_home.append("g")
                      .attr("transform", "translate(0,-15)");

    first_legend_svg.append("text")
              .attr('y', 5)
              .attr('x', 40)
              .style("font", "14px sans-serif")
              .style("font-weight", "bold");

    first_legend_svg.append("text")
              .attr('x', 30)
              .attr('y', 33)
              .style("font", "14px sans-serif")
              .style("font-weight", "bold")
              .text(startValue);

    first_legend_svg.append('rect')
              .attr('x', 40)
              .attr('y', 20) // Adjust the y-coordinate to add space between gradient and text box
              .attr('width', (width*0.9))
              .attr('height', 20)
              .style('fill', 'url(#grad_' + score + ')');

    first_legend_svg.append("text")
              .attr('x', width + 20)
              .attr('y', 33)
              .style("font", "14px sans-serif")
              .style("font-weight", "bold")
              .text(endValue);

    legendLabel = first_legend_svg.select('text')['_groups'][0][0].getBBox().width*1.05 + 0.5 * 10;

    first_legend_svg.select("text")
              .attr("x", (width + margin.left + margin.right - legendLabel)/2);

    // Draw the black vertical line at the defined value
    var lineX = (definedValue - startValue) / (endValue - startValue) * (width * 0.9) + 40;
    svg_home.append("line")
            .attr("x1", lineX)
            .attr("y1", height - 45) // Adjust the y-coordinate of the top arrowhead 90
            .attr("x2", lineX)
            .attr("y2", height - 25) // Adjust the y-coordinate of the bottom arrowhead 125
            .attr("stroke", "black")
            .attr("stroke-width", 2);

    // Add arrowheads
    svg_home.append("path")
            .attr("d", "M " + lineX + " " + (height - 45) + " L " + (lineX - 5) + " " + (height - 50) + " L " + (lineX + 5) + " " + (height - 50) + " Z")
            .attr("fill", "black");

    svg_home.append("path")
            .attr("d", "M " + lineX + " " + (height - 25) + " L " + (lineX - 5) + " " + (height - 20) + " L " + (lineX + 5) + " " + (height - 20) + " Z")
            .attr("fill", "black");

}

// This function generates two concentric circles made of aminoAcids
// then creates lines connecting inner and outer aminoacids, given interactions between them

function draw_interactions_in_circles(location, interactions, inner_data, outer_data, conversion_dict, signalling_protein) {
  // D3 select the SVG
  const svg2 = d3.select("#" + location);

  // Types of interactions and their colors
  const interactionTypes = {
    'Aromatic': '#689235',
    'Hydrophobic': '#A6E15F',
    'Ionic': '#005693',
    'Polar': '#7030a0',
    'Van der waals': '#d9d9d9'
  };

  // Continuous Line: ""
  // Short Dashes: "4,2"
  // Long Dashes: "8,2"
  // Dots: "1,3"
  // Dash-Dot-Dash: "8,2,1,2"

  const strokeShape = {
    'Aromatic': '',
    'Hydrophobic': '4,2',
    'Ionic': '8,2',
    'Polar': '1,3',
    'Van der waals': '8,2,1,2'
    }

  const allColors = {
    'Aromatic': '#689235',
    'Hydrophobic': '#A6E15F',
    'Ionic': '#005693',
    'Polar': '#7030a0',
    'Van der waals': '#d9d9d9',
    'TM1': '#000000',  //Black
    'ICL1': '#0E0E0E',  //Very Dark Gray
    'TM2': '#1C1C1C',  //Darker Gray
    'ICL2': '#2A2A2A',  //Dark Gray
    'TM3': '#383838',  //Medium-Dark Gray
    'ICL3': '#464646',  //Medium Gray
    'TM4': '#545454',  //Medium Gray
    'TM5': '#626262',  //Medium-Light Gray
    'TM6': '#707070',  //Lighter Gray
    'TM7': '#7E7E7E',  //Light Gray
    'ICL4': '#8C8C8C',  //Light Gray
    'H8': '#9A9A9A',   //Very Light Gray
    'C-term': '#CCCCCC',  //Light Gray
    'G.HN': '#FF0000',
    'G.hns1': '#FF0C00',
    'G.S1': '#FF1900',
    'G.s1h1': '#FF2500',
    'G.H1': '#FF3200',
    'G.h1ha': '#FF3E00',
    'H.HA': '#FF4B00',
    'H.hahb': '#FF5700',
    'H.HB': '#FF6400',
    'H.hbhc': '#FF7100',
    'H.HC': '#FF7D00',
    'H.hchd': '#FF8A00',
    'H.HD': '#FF9600',
    'H.hdhe': '#FFA300',
    'H.HE': '#FFAF00',
    'H.hehf': '#FFBC00',
    'H.HF': '#FFC900',
    'G.hfs2': '#FFD500',
    'G.S2': '#FFE200',
    'G.s2s3': '#FFEE00',
    'G.S3': '#FFFB00',
    'G.s3h2': '#F9F300',
    'G.H2': '#F4EC00',
    'G.h2s4': '#F0E400',
    'G.S4': '#EBDD00',
    'G.s4h3': '#E7D500',
    'G.H3': '#E2CE00',
    'G.h3s5': '#DEC600',
    'G.S5': '#D9BF00',
    'G.s5hg': '#D5B700',
    'G.HG': '#D0B000',
    'G.hgh4': '#CCA800',
    'G.H4': '#C7A100',
    'G.h4s6': '#C39900',
    'G.S6': '#BE9200',
    'G.s6h5': '#BA8A00',
    'G.H5': '#876400',
    "ns1": "#006400",
    "S1": "#006A00",
    "s1s2": "#007100",
    "S2": "#007700",
    "s2s3": "#007D00",
    "S3": "#008400",
    "s3s4": "#008A00",
    "S4": "#009100",
    "s4s5": "#009700",
    "S5": "#009D00",
    "s5s6": "#00A400",
    "S6": "#00AA00",
    "s6h1": "#00B100",
    "H1": "#00B700",
    "h1s7": "#00BD00",
    "S7": "#00C400",
    "s7s8": "#00CA00",
    "S8": "#00D100",
    "s8s9": "#00D700",
    "S9": "#00DD00",
    "s9s10": "#00E400",
    "S10": "#00EA00",
    "s10s11": "#00F100",
    "S11": "#00F700",
    "s11s12": "#00FD00",
    "S12": "#00FF00",
    "s12s13": "#00FF06",
    "S13": "#00FF0D",
    "s13s14": "#00FF13",
    "S14": "#00FF1A",
    "s14s15": "#00FF20",
    "S15": "#00FF27",
    "s15s16": "#00FF2D",
    "S16": "#00FF34",
    "s16s17": "#00FF3A",
    "S17": "#00FF40",
    "s17s18": "#00FF47",
    "S18": "#00FF4D",
    "s18s19": "#00FF54",
    "S19": "#00FF5A",
    "s19s20": "#00FF61",
    "S20": "#00FF67",
    "s20c": "#00FF6D"
  };

  const aminoAcids = [
    { aminoAcid: "Alanine", singleLetter: "A" },
    { aminoAcid: "Arginine", singleLetter: "R" },
    { aminoAcid: "Asparagine", singleLetter: "N" },
    { aminoAcid: "Aspartic acid", singleLetter: "D" },
    { aminoAcid: "Cysteine", singleLetter: "C" },
    { aminoAcid: "Glutamine", singleLetter: "Q" },
    { aminoAcid: "Glutamic acid", singleLetter: "E" },
    { aminoAcid: "Glycine", singleLetter: "G" },
    { aminoAcid: "Histidine", singleLetter: "H" },
    { aminoAcid: "Isoleucine", singleLetter: "I" },
    { aminoAcid: "Leucine", singleLetter: "L" },
    { aminoAcid: "Lysine", singleLetter: "K" },
    { aminoAcid: "Methionine", singleLetter: "M" },
    { aminoAcid: "Phenylalanine", singleLetter: "F" },
    { aminoAcid: "Proline", singleLetter: "P" },
    { aminoAcid: "Serine", singleLetter: "S" },
    { aminoAcid: "Threonine", singleLetter: "T" },
    { aminoAcid: "Tryptophan", singleLetter: "W" },
    { aminoAcid: "Tyrosine", singleLetter: "Y" },
    { aminoAcid: "Valine", singleLetter: "V" }
  ];

  const presetColors = {'D': ['#E60A0A', '#FDFF7B'],'E': ['#E60A0A', '#FDFF7B'],
                        'K': ['#145AFF', '#FDFF7B'],'R': ['#145AFF', '#FDFF7B'],
                        'S': ['#A70CC6', '#FDFF7B'],'T': ['#A70CC6', '#FDFF7B'],
                        'N': ['#A70CC6', '#FDFF7B'],'Q': ['#A70CC6', '#FDFF7B'],
                        'V': ['#E6E600', '#000000'],'L': ['#E6E600', '#000000'],
                        'I': ['#E6E600', '#000000'],'A': ['#E6E600', '#000000'],
                        'M': ['#E6E600', '#000000'],'F': ['#18FF0B', '#000000'],
                        'Y': ['#18FF0B', '#000000'],'W': ['#0BCF00', '#000000'],
                        'H': ['#0093DD', '#000000'],'P': ['#CC0099', '#FDFF7B'],
                        'C': ['#B2B548', '#000000'],'G': ['#FF00F2', '#000000'],
                        '-': ['#FFFFFF', '#000000'],'+': ['#FFFFFF', '#000000']};


  // Extract all 'segment' values from the data array
  const innerArray = inner_data.map(item => item.segment);
  const outerArray = outer_data.map(item => item.segment);

  // Remove duplicates by converting to a Set and then back to an array
  const inner_legend = Array.from(new Set(innerArray));
  const outer_legend = Array.from(new Set(outerArray));

  const colorOrder = Object.keys(allColors);
  outer_legend.sort((a, b) => colorOrder.indexOf(a) - colorOrder.indexOf(b));
  inner_legend.sort((a, b) => colorOrder.indexOf(a) - colorOrder.indexOf(b));

  const countInner = inner_data.length;
  const countOuter = outer_data.length;
  // Generate arrays for inner and outer circles
  const beadInfo = {
    innerCircle: inner_data,
    outerCircle: outer_data
  };

  const numInnerBeads = beadInfo.innerCircle.length;

  const bead_gap = beadInfo.outerCircle.length - beadInfo.innerCircle.length;

  function sanitizeClassName(name) {
    return name.toLowerCase().replace(/\s+/g, '-');
  }

  const tooltip = d3.select("body").append("div")
    .attr("class", "tooltip")
    .style("opacity", 0);

  // Function to create a bead
  function createBead(cx, cy, aminoAcid, segment, index, circleType, beadRadius, centroidX, centroidY, grn_number, list_int, seq_number) {
    // Tooltip show function
    function showTooltip(d) {
      const event = window.event ? window.event : d3.event;
      tooltip.transition()
        .duration(200)
        .style("opacity", .9);
      tooltip.html("Amino Acid: " + aminoAcid + "<br/>" +
                  "Segment: " + segment + "<br/>" +
                  "Generic residue number: " + grn_number + "<br/>" +
                  "Sequence number: " + seq_number + "<br/>" +
                  "Interactions: " + list_int + "<br/>"
                )
        .style("left", (event.pageX + 5) + "px")
        .style("top", (event.pageY - 28) + "px");
    }

    // Tooltip hide function
    function hideTooltip(d) {
      tooltip.transition()
        .duration(500)
        .style("opacity", 0);
    }

    let strokeWidth;

    if (circleType == 'outer'){
      strokeWidth = 5;
    } else {
      strokeWidth = 3;
    };

    const bead = svg2.append("circle")
      .attr("cx", cx)
      .attr("cy", cy)
      .attr("r", beadRadius)
      .attr("fill", "white")
      .attr("stroke", allColors[segment])
      .attr("stroke-width", strokeWidth)
      .attr("data-index", index)
      .attr("data-circle", circleType)
      .attr("data-aa", aminoAcid)
      .attr("data-segment", segment)
      .classed(segment, true)
      .attr("id", `bead-${circleType}-${index}`)
      .on("click", function(d) {
        // Clear all previous highlights
        d3.selectAll("path").filter(function() {
          return !d3.select(this).attr("segment");
        })
        .each(function(d, i) {
          d3.select(this).attr("stroke-opacity", 0);
        });
        d3.selectAll(`path[data-${circleType}-bead="${index}"]`).attr("stroke-opacity", 1);
      })
      .on("mouseover", showTooltip)
      .on("mouseout", hideTooltip);

    // Add text inside the circle to represent the amino acid
    svg2.append("text")
      .attr("x", cx)
      .attr("y", cy)
      .attr("dy", 9)
      .attr("text-anchor", "middle")
      .text(aminoAcid)
      .style("font-size", "23px")
      .on("mouseover", showTooltip)
      .on("mouseout", hideTooltip);

    // Calculate radial text position
    const angleToCentroid = Math.atan2(cy - centroidY, cx - centroidX);

    const textDistanceFromCenter = beadRadius * 2 + 5; // distance from the center of the bead

    let textX = cx + textDistanceFromCenter * Math.cos(angleToCentroid * (Math.PI / 180));
    let textY = cy + textDistanceFromCenter * Math.sin(angleToCentroid * (Math.PI / 180));

    // Determine text rotation and alignment
    let textRotation = angleToCentroid * (180 / Math.PI)
    let textAnchor = "start";

    // Correct the orientation for better readability
    if (circleType == 'inner') {
      inner_number = grn_number.split('.')[2];
      if (textRotation > 90 || textRotation < -90 ) {
        textRotation += 180;
        textAnchor = "end";
      } else{
        textAnchor = "start";
        textX = cx + (textDistanceFromCenter - 4 * beadRadius - 10) * Math.cos(angleToCentroid * (Math.PI / 180));
        textY = cy + (textDistanceFromCenter - 2 * beadRadius) * Math.sin(angleToCentroid * (Math.PI / 180));
      }
    } else {
        inner_number= grn_number.replace(/\..*x/, 'x');
        // Correct the orientation for better readability
        if (textRotation > 90 || textRotation < -90 ) {
          textRotation += 180;
          textAnchor = "end";  // Switch to "end"

          // Shift the text position so it starts from the opposite side
          textX = cx + (textDistanceFromCenter - 4 * beadRadius + 10) * Math.cos(angleToCentroid * (Math.PI / 180));
          textY = cy + (textDistanceFromCenter - 2 * beadRadius) * Math.sin(angleToCentroid * (Math.PI / 180));
        } else {
          textX = textX - 20;
        }
    }

    let generic_label, protein_label;
    if(circleType=='inner'){
      generic_label = grn_number + "_gprot";
      protein_label = seq_number + "_gprot";
    } else {
      generic_label = grn_number + "_GPCR";
      protein_label = seq_number + "_GPCR";
    }
    // Create text
    svg2.append("text")
        .attr("x", textX)
        .attr("y", textY)
        .attr("dy", "0.3em")  // Adjust as needed for vertical alignment
        .attr("text-anchor", textAnchor)
        .attr("transform", `rotate(${textRotation},${cx},${cy})`) // .attr("transform", `rotate(${angleToCentroid * (180 / Math.PI)},${cx},${cy})`)
        .text(inner_number)
        .attr("generic-nr", generic_label)
        .attr("protein-nr", protein_label)
        .style("font-size", "18px")
        .style("font-family", "Palatino")
        .attr("font-weight", "bold")
        .style("fill", "#111");  // or any color of your choice
  }

  // Function to create a circle of beads
  function createCircle(cx, cy, beadRadius, beads, circleType, circleRadius) {
    const numBeads = beads.length;

    // Add additional distance for the labels, assuming each label needs 15 units of space
    const internalCircleRadius = circleRadius - 35;

    // The centroid of the circle is simply its center, denoted by (cx, cy)
    const centerX = cx;
    const centerY = cy;

    let lastSegment = null;
    let textAnchor;
    let segmentStartAngle = 0;
    const gapInRadians = 10 / internalCircleRadius; // 10 pixels space

    beads.forEach((bead, index) => {
      const angle = (index / numBeads) * 2 * Math.PI;
      const x = cx + circleRadius * Math.cos(angle);
      const y = cy + circleRadius * Math.sin(angle);

      if (circleType == 'inner') {
        if (lastSegment !== bead.segment) {
          if (lastSegment !== null) {
            const adjustedStartAngle = segmentStartAngle + gapInRadians;
            const adjustedEndAngle = angle - gapInRadians;

            // Draw the segment label at the midpoint of the arc segment
            let midAngle = (adjustedStartAngle + adjustedEndAngle) / 2;
            let labelX = cx + (internalCircleRadius - 55) * Math.cos(midAngle);
            let labelY = cy + (internalCircleRadius - 55) * Math.sin(midAngle);

            let textAngleDeg = (midAngle * 180 / Math.PI);  // Convert to degrees and align the text

            if (textAngleDeg > 90 || textAngleDeg < -90 ) {
              textAngleDeg += 180;
              textAnchor = "end";
            } else{
              textAnchor = "start";
              // labelX = cx + (internalCircleRadius - 4 * beadRadius - 10) * Math.cos(midAngle * (Math.PI / 180));
              // labelY = cy + (internalCircleRadius - 2 * beadRadius) * Math.sin(midAngle * (Math.PI / 180));
            }

            svg2.append("text")
                .attr("x", labelX)
                .attr("y", labelY)
                .attr("text-anchor", textAnchor)
                // .attr("alignment-baseline", "middle")
                .attr("transform", `rotate(${textAngleDeg}, ${labelX}, ${labelY})`)
                .text(lastSegment)
                .style("font-size", "14px")
                .attr("font-weight", "bold")
                .style("font-family", "Palatino")
                .style("fill", "#333");

            // Draw the external circle segment
            // drawExternalSegment(cx, cy, internalCircleRadius, adjustedStartAngle, adjustedEndAngle, lastSegment);
          }

          segmentStartAngle = angle;
          lastSegment = bead.segment;
        }
      }

      createBead(x, y, bead.aminoAcid, bead.segment, index, circleType, beadRadius, centerX, centerY, bead.generic_number, bead.interactions, bead.sequence_number);
    });

    // Draw the remaining external circle segment here, if needed
    const adjustedStartAngle = segmentStartAngle + gapInRadians;
    const adjustedEndAngle = 2 * Math.PI - gapInRadians;
    let smallTextAnchor;
    // Draw the segment label for the remaining segment
    const midAngle = (adjustedStartAngle + adjustedEndAngle) / 2;
    let labelX = cx + (internalCircleRadius - 45) * Math.cos(midAngle);
    let labelY = cy + (internalCircleRadius - 45) * Math.sin(midAngle);

    let textAngleDeg = (midAngle * 180 / Math.PI);  // Convert to degrees and align the text

    if (textAngleDeg > 90 || textAngleDeg < -90 ) {
      textAngleDeg += 180;
      smallTextAnchor = "end";
    } else{
      smallTextAnchor = "start";
    }

    svg2.append("text")
        .attr("x", labelX)
        .attr("y", labelY)
        .attr("text-anchor", smallTextAnchor)
        // .attr("alignment-baseline", "middle")
        .attr("transform", `rotate(${textAngleDeg}, ${labelX}, ${labelY})`)
        .text(lastSegment)
        .style("font-size", "14px")
        .attr("font-weight", "bold")
        .style("font-family", "Palatino")
        .style("fill", "#333");

    let lines;
    let segment_label;
    // Split the centerText into two lines
    if (signalling_protein == 'af-signprot') {
      lines = 'G\nprotein'.split('\n');
      segment_label = 'G Protein segments:';
    } else {
      lines = 'Arrestin'.split('\n');
      segment_label = 'Arrestin segments:';
    };


    // Determine the vertical offset for the text lines to make them appear centered
    const fontSize = 14; // Assuming a font size of 20 for illustration, adjust as needed
    const lineHeight = fontSize * 1.2;  // Assuming a line height 20% greater than font size, adjust as needed
    const totalHeight = lines.length * lineHeight;

    // Draw text lines in the center of the circle
    let maxWidth = 0;
    let maxHeight = 0;

    lines.forEach((line, index) => {
        const textElement = svg2.append("text")
            .attr("x", cx)  // set x-coordinate to the center x-coordinate
            .attr("y", cy - totalHeight / 2 + (index + 0.5) * lineHeight) // adjust y-coordinate based on line index
            .attr("text-anchor", "middle")  // align horizontally to the middle
            .attr("alignment-baseline", "middle")  // align vertically to the middle
            .text(line)  // set the text string for the current line
            .style("font-size", `${fontSize}px`)  // set font size
            .attr("font-weight", "bold")  // set font weight
            .style("font-family", "Arial")  // set font family
            .style("fill", "white");  // set font color

        // Create a bounding box for the text
        const textBoundingBox = textElement.node().getBBox();

        // Update maxWidth and maxHeight based on the current bounding box
        maxWidth = Math.max(maxWidth, textBoundingBox.width);
        maxHeight = Math.max(maxHeight, textBoundingBox.height);
    });

    // Create a single rectangle that serves as the background for all text lines
    svg2.insert("rect", ":first-child")
        .attr("x", cx - maxWidth / 2 - 5) // Adjust the padding as needed
        .attr("y", cy - totalHeight/2 - maxHeight/2 +5) // Adjust the padding as needed
        .attr("width", maxWidth + 10) // Use the maximum width
        .attr("height", totalHeight + 4) // Use the total height of all text lines
        .attr("rx", 6) // Adjust the border radius as needed
        .attr("ry", 6) // Adjust the border radius as needed
        .style("fill", "gray"); // set background color to gray
  }

  function averageAngle(angles) {
      let sumX = 0;
      let sumY = 0;

      for (let angle of angles) {
          sumX += Math.cos(angle);
          sumY += Math.sin(angle);
      }

      return Math.atan2(sumY / angles.length, sumX / angles.length);
  }

  // Function to create outer beads based on the connections dictionary
  function createOuterBeads(centroidX, centroidY, innerCircleRadius, outerCircleRadius, beadConnections, beads) {
      const beadRadius = 20; // Or whatever the radius of the outer beads is
      const placedBeads = []; // To store the angles at which beads have been placed

      for (const [outerBeadId, connections] of Object.entries(beadConnections)) {
          const connectionAngles = connections.map(connId => {
              const innerBead = document.querySelector(`#bead-inner-${connId}`);
              const beadCX = parseFloat(innerBead.getAttribute('cx'));
              const beadCY = parseFloat(innerBead.getAttribute('cy'));
              return Math.atan2(beadCY - centroidY, beadCX - centroidX);
          });

          connectionAngles.sort((a, b) => a - b);

          // Find the maximum angle difference
          let maxDifference = -Infinity;
          let chosenAnglePair = [];
          for (let i = 0; i < connectionAngles.length; i++) {
              let diff = connectionAngles[(i+1) % connectionAngles.length] - connectionAngles[i];
              if (diff < 0) {
                  diff += 2 * Math.PI;
              }
              if (diff > maxDifference) {
                  maxDifference = diff;
                  chosenAnglePair = [connectionAngles[i], connectionAngles[(i+1) % connectionAngles.length]];
              }
          }
          // Calculate average vector direction
          const avgCos = connectionAngles.reduce((sum, angle) => sum + Math.cos(angle), 0) / connectionAngles.length;
          const avgSin = connectionAngles.reduce((sum, angle) => sum + Math.sin(angle), 0) / connectionAngles.length;

          // let avgAngle = (chosenAnglePair[0] + chosenAnglePair[1]) / 2;
          // let avgAngle = connectionAngles.reduce((a, b) => a + b, 0) / connectionAngles.length;
          // let avgAngle = Math.atan2(avgSin, avgCos);
          let avgAngle = averageAngle(chosenAnglePair);

          // Correction factor for overlapping beads
          let angleCorrection = 0.08;  // This is the angular shift applied. Adjust as needed.
          while (placedBeads.some(angle => Math.abs(angle - avgAngle) < angleCorrection)) {
              avgAngle += angleCorrection;
          }
          placedBeads.push(avgAngle);

          const outerBeadX = centroidX + outerCircleRadius * Math.cos(avgAngle);
          const outerBeadY = centroidY + outerCircleRadius * Math.sin(avgAngle);

          let bead = beads[outerBeadId];
          // You can use the createBead function here
          createBead(outerBeadX, outerBeadY, bead.aminoAcid, bead.segment, outerBeadId, 'outer', beadRadius, centroidX, centroidY, bead.generic_number, bead.interactions, bead.sequence_number);
      }
  };

  function computeBezierPath(startX, startY, endX, endY, beadRadius) {
      // Calculate the direction vector between start and end
      let dx = endX - startX;
      let dy = endY - startY;

      // Calculate the normalized direction vector
      let len = Math.sqrt(dx * dx + dy * dy);
      let nx = dx / len;
      let ny = dy / len;

      // Compute the points where the curve touches the bead border
      let adjustedStartX = startX + nx * beadRadius;
      let adjustedStartY = startY + ny * beadRadius;
      let adjustedEndX = endX - nx * beadRadius;
      let adjustedEndY = endY - ny * beadRadius;

      // Calculate the normal vector (perpendicular)
      let nxPerpendicular = -ny;
      let nyPerpendicular = nx;

      // Determine a control point distance factor for the curve
      // We'll give a gentle pull for the control point near the destination
      let curveFactor = 0.1 * len;

      // Determine the control points
      // The first control point will be in the direct line between start and end
      let controlPoint1X = adjustedStartX + 0.5 * (adjustedEndX - adjustedStartX);
      let controlPoint1Y = adjustedStartY + 0.5 * (adjustedEndY - adjustedStartY);
      // The second control point will be slightly pulled out using the normal vector
      let controlPoint2X = adjustedEndX - nx * curveFactor + nxPerpendicular * curveFactor;
      let controlPoint2Y = adjustedEndY - ny * curveFactor + nyPerpendicular * curveFactor;

      // Create the Bezier curve path
      let pathData = `M${adjustedStartX},${adjustedStartY} C${controlPoint1X},${controlPoint1Y} ${controlPoint2X},${controlPoint2Y} ${adjustedEndX},${adjustedEndY}`;

      return pathData;
  }

  function drawConnectionLines(connections, svg2, innerCircleRadius) {
      connections.forEach(connection => {
          const innerBead = d3.select(`#bead-inner-${connection.innerIndex}`);
          const outerBead = d3.select(`#bead-outer-${connection.outerIndex}`);

          const innerBeadX = parseFloat(innerBead.attr('cx'));
          const innerBeadY = parseFloat(innerBead.attr('cy'));
          const outerBeadX = parseFloat(outerBead.attr('cx'));
          const outerBeadY = parseFloat(outerBead.attr('cy'));

          const beadRadius = 20;
          const pathData = computeBezierPath(outerBeadX, outerBeadY, innerBeadX, innerBeadY, beadRadius);

          const sanitizedType = sanitizeClassName(connection.type);

          svg2.append("path")
              .attr("d", pathData)
              .attr("id", `path-${connection.innerIndex}-${connection.outerIndex}`)
              .attr("data-inner-bead", connection.innerIndex)
              .attr("data-outer-bead", connection.outerIndex)
              .attr("stroke", "black")
              .attr("stroke-width", 2)
              .attr("fill", "none")
              .attr("interaction", connection.type)
              .attr('inner-chain', connection.innerChain)
              .attr('outer-chain', connection.outerChain)
              .attr("stroke-dasharray", strokeShape[connection.type])
              .classed(sanitizedType, true);
      });
  }


  // Function to draw external circle segment
  function drawExternalSegment(cx, cy, radius, startAngle, endAngle, segment) {
    const startX = cx + radius * Math.cos(startAngle);
    const startY = cy + radius * Math.sin(startAngle);
    const endX = cx + radius * Math.cos(endAngle);
    const endY = cy + radius * Math.sin(endAngle);

    const d = [
      `M ${startX} ${startY}`,  // Move to the start point
      `A ${radius} ${radius} 0 0 1 ${endX} ${endY}`  // Draw an arc to the end point
    ].join(" ");

    svg2.append("path")
      .attr("d", d)
      .attr("stroke", allColors[segment])
      .attr("interaction", segment)
      .attr("segment", "external")
      .attr("stroke-width", 4)
      .attr("fill", "none");
  }

  function calculateCentroids(d3Selection) {
    let sumX = 0;
    let sumY = 0;
    let count = 0;

    d3Selection.each(function() {
      const circle = d3.select(this);
      sumX += +circle.attr('cx');
      sumY += +circle.attr('cy');
      count++;
    });

    const centroidX = sumX / count;
    const centroidY = sumY / count;

    return { centroidX, centroidY };
  }

  function createConnection(innerIndex, outerIndex, type, countInner, beadRadius, centroidX, centroidY, innerChain, outerChain, smallRadius, bigRadius) {
    const innerBead = d3.select(`#bead-inner-${innerIndex}`);
    const outerBead = d3.select(`#bead-outer-${outerIndex}`);

    // Calculate direction vector for inner bead
    const innerX = +innerBead.attr("cx");
    const innerY = +innerBead.attr("cy");
    const innerDirX = innerX - centroidX;
    const innerDirY = innerY - centroidY;

    // Calculate direction vector for outer bead
    const outerX = +outerBead.attr("cx");
    const outerY = +outerBead.attr("cy");
    const outerDirX = outerX - centroidX;
    const outerDirY = outerY - centroidY;

    const lengthInner = Math.sqrt(innerDirX * innerDirX + innerDirY * innerDirY);
    const unitInnerX = -innerDirX / lengthInner;
    const unitInnerY = -innerDirY / lengthInner;

    const lengthOuter = Math.sqrt(outerDirX * outerDirX + outerDirY * outerDirY);
    const unitOuterX = -outerDirX / lengthOuter;
    const unitOuterY = -outerDirY / lengthOuter;

    // Calculate the start and end points of the bead circles
    // For the start point: on the external border of the outer bead, facing INNER_BORDER_RADIUS
    const startX = outerX + unitOuterX * beadRadius;
    const startY = outerY + unitOuterY * beadRadius;

    // For the end point: on the external border of the inner bead, facing INNER_BORDER_RADIUS
    const endX = innerX - unitInnerX * beadRadius;
    const endY = innerY - unitInnerY * beadRadius;

    const angleToCentroidInner = Math.atan2(innerY - centroidY, innerX - centroidX);
    const angleToCentroidOuter = Math.atan2(outerY - centroidY, outerX - centroidX);

    const INNER_BORDER_RADIUS = smallRadius + 22;
    const MIDDLE_BORDER_RADIUS = smallRadius + 3 * beadRadius;
    const OUTER_BORDER_RADIUS = bigRadius - 4 * beadRadius;

    const scaleFactor = 1.4;

    const tangentialOuterX = centroidX + Math.cos(angleToCentroidOuter) * OUTER_BORDER_RADIUS;
    const tangentialOuterY = centroidY + Math.sin(angleToCentroidOuter) * OUTER_BORDER_RADIUS;

    const dynamicScaleFactor = (countInner > 15) ? 1.6 : 1.4; // Adjust based on your requirements
    let tangentialMiddleX = centroidX + Math.cos((angleToCentroidInner + angleToCentroidOuter) / 2) * MIDDLE_BORDER_RADIUS * dynamicScaleFactor;
    let tangentialMiddleY = centroidY + Math.sin((angleToCentroidInner + angleToCentroidOuter) / 2) * MIDDLE_BORDER_RADIUS * dynamicScaleFactor;

    let tangentialInnerX = centroidX + Math.cos(angleToCentroidInner) * INNER_BORDER_RADIUS;
    let tangentialInnerY = centroidY + Math.sin(angleToCentroidInner) * INNER_BORDER_RADIUS;

    // Function to check if a given point overlaps with any inner bead
    function overlapsWithInnerBead(x, y) {
        for (let i = 0; i < countInner; i++) {
            if (i !== innerIndex) {
                const otherBead = d3.select(`#bead-inner-${i}`);
                const otherX = +otherBead.attr("cx");
                const otherY = +otherBead.attr("cy");
                const distance = Math.sqrt((x - otherX) ** 2 + (y - otherY) ** 2);
                if (distance < 2 * beadRadius) {
                    return true;
                }
            }
        }
        return false;
    }

    let adjustedEndAngle = Math.atan2(tangentialInnerY - innerY, tangentialInnerX - innerX);
    let adjustedEndX = innerX + Math.cos(adjustedEndAngle) * beadRadius;
    let adjustedEndY = innerY + Math.sin(adjustedEndAngle) * beadRadius;

    if (overlapsWithInnerBead(adjustedEndX, adjustedEndY)) {
        const offsetAngle = 0.05 * innerIndex;
        adjustedEndAngle += offsetAngle;
        adjustedEndX = innerX + Math.cos(adjustedEndAngle) * beadRadius;
        adjustedEndY = innerY + Math.sin(adjustedEndAngle) * beadRadius;
    }

    let pathData;

    const isClockwise = angleToCentroidOuter > angleToCentroidInner;

    let controlPoint1X = startX + (tangentialOuterX - startX) / 2;
    let controlPoint1Y = startY + (tangentialOuterY - startY) / 2;
    let controlPoint2X = endX + (tangentialInnerX - endX) / 2;
    let controlPoint2Y = endY + (tangentialInnerY - endY) * 20;

    if(countInner > 15){
      if(isClockwise){
        controlPoint1X = startX + (tangentialMiddleX - startX) / 2;
        controlPoint1Y = startY + (tangentialMiddleY - startY) / 2;
        controlPoint2X = endX + (tangentialInnerX - endX) / 2;
        controlPoint2Y = endY + (tangentialInnerY - endY) / 2;
      } else {
        controlPoint1X = (startX + tangentialMiddleX) / 2;
        controlPoint1Y = (startY + tangentialMiddleY) / 2;
        controlPoint2X = (endX + tangentialInnerX) / 2;
        controlPoint2Y = (endY + tangentialInnerY) / 2;
      }
    }

    if(countInner < 15){
      pathData = [
        `M ${startX} ${startY}`,
        `C ${controlPoint1X} ${controlPoint1Y}, ${controlPoint2X} ${controlPoint2Y}, ${adjustedEndX} ${adjustedEndY}`
      ].join(" ");
    } else {
      pathData = [
        `M ${startX} ${startY}`,
        `C ${controlPoint1X} ${controlPoint1Y}, ${tangentialMiddleX} ${tangentialMiddleY}, ${tangentialInnerX} ${tangentialInnerY}`,
        `C ${controlPoint2X} ${controlPoint2Y}, ${controlPoint2X} ${controlPoint2Y}, ${adjustedEndX} ${adjustedEndY}`
      ].join(" ");
    }

    const sanitizedType = sanitizeClassName(type);

    // Draw the path
    svg2.append("path")
      .attr("d", pathData)
      .attr("id", `path-${innerIndex}-${outerIndex}`)
      .attr("data-inner-bead", innerIndex)
      .attr("data-outer-bead", outerIndex)
      .attr("stroke", "black")
      .attr("stroke-width", 2)  // Border width
      .attr("fill", "none")
      .attr("interaction", type)
      .attr('inner-chain', innerChain)
      .attr('outer-chain', outerChain)
      .attr("stroke-dasharray", strokeShape[type])
      .classed(sanitizedType, true);
  }

  let smallRadius;
  let bigRadius;

  if (countInner <= 15){
    smallRadius = 25 * 20 / Math.PI * 0.85;
  } else if (countInner > 15 && countInner <= 35) {
    smallRadius = 35 * 20 / Math.PI * 0.85;
  } else {
    smallRadius = countInner * 20 / Math.PI * 0.85;
  };

  if (countInner <= 15){
    bigRadius = smallRadius * 1.8;
  } else if (countInner > 15 && countInner <= 35 ){
    bigRadius = smallRadius * 2.2;
  } else {
    bigRadius = (countOuter * 20 / Math.PI * 0.95) * 2.2;
  }

  let svgH = (bigRadius * 2) + 400;
  svg2.attr("height", svgH);
  svg2.attr("width", svgH);

  // For example, with a bead radius of 20
  createCircle(bigRadius*1.3, bigRadius*1.3, 20, beadInfo.innerCircle, 'inner', smallRadius);
  createOuterBeads(bigRadius*1.3, bigRadius*1.3, smallRadius, bigRadius, conversion_dict, beadInfo.outerCircle);

  drawConnectionLines(interactions, svg2, smallRadius);
  const innerBeadSelection = svg2.selectAll("circle[data-circle='inner']");  // Assuming the inner beads have a class 'inner-bead'
  const { centroidX, centroidY } = calculateCentroids(innerBeadSelection);

  // interactions.forEach(({ innerIndex, outerIndex, type, innerChain, outerChain }) => createConnection(innerIndex, outerIndex, type, countInner, 20, centroidX, centroidY, innerChain, outerChain, smallRadius, bigRadius));


  // Generate legend
  const outerCircleLabel = svg2.append("g")
    .attr("transform", `translate(0, ${bigRadius*1.3})`);

  outerCircleLabel.append("rect")  // Add a rectangle for the background
      .attr("x", 1)
      .attr("y", -14) // Adjust the position as needed
      .attr("width", 50) // Adjust the width as needed
      .attr("height", 25) // Adjust the height as needed
      .attr("rx", 6) // Adjust the border radius as needed
      .attr("ry", 6) // Adjust the border radius as needed
      .attr("fill", "gray"); // Set the background color to gray

  outerCircleLabel.append("text")
                  .attr("x", 25)
                  .attr("y", 0)
                  .text("GPCR")
                  .attr("text-anchor", "middle")  // align horizontally to the middle
                  .attr("alignment-baseline", "middle")  // align vertically to the middle
                  .attr("font-family", "Arial")
                  .attr("font-size", "14px")
                  .attr("font-weight", "bold")
                  .attr("fill", "white");

  // Generate legend
  const interactionLegend = svg2.append("g")
    .attr("transform", `translate(72, ${svgH - 60})`);

  // Generate inner legend
  const innercircleLegend = svg2.append("g")
    .attr("transform", `translate(72, ${svgH - 30})`);

  // Generate outer legend
  const outercircleLegend = svg2.append("g")
    .attr("transform", `translate(72, ${svgH})`);

  // Header for the interaction legend
  const headerInteractions = interactionLegend.append("text")
                              .attr("x", 0)
                              .attr("y", -20)
                              .attr("class", "legend-header")
                              .text("Interactions:")
                              .attr("font-family", "Arial")
                              .attr("font-size", "16px")
                              .attr("font-weight", "bold")
                              .attr("fill", "black");

  const txtLen = headerInteractions.node().getComputedTextLength();

  let cumulativeX = txtLen + 20;

  Object.entries(interactionTypes).forEach(([type, color], i) => {
    const legendItem = interactionLegend.append("g")
      .attr("transform", `translate(${cumulativeX}, 0)`);

    // Draw the rectangle at y=0
    legendItem.append("rect")
      .attr("width", 20)
      .attr("height", 20)
      .attr("y", -35)  // <-- Set y to 0 for the rectangle
      .attr("fill", color)
      .attr("stroke", "black")
      .attr("stroke-width", 1)
      .on("click", function() {
        if(d3.event) {
          d3.event.stopPropagation();
        }
        const sanitizedType = sanitizeClassName(type);
        d3.selectAll("path").attr("stroke-opacity", 0);
        d3.selectAll("path[segment='external']").attr("stroke-opacity", 1);
        d3.selectAll(`.${sanitizedType}`).attr("stroke-opacity", 1);
      });

    // Draw the text at y=15
    const textElement = legendItem.append("text")
      .attr("x", 25)
      .attr("y", -20)  // <-- Set y to 15 for the text
      .text(type);

    // Get text width and update the cumulativeX for the next legend item
    const textWidth = textElement.node().getComputedTextLength();
    cumulativeX += textWidth + 50; // 50 is the space between the rectangle and the next item
  });

  // Creating the inner legend structur

  let segment_label;
  // Split the centerText into two lines
  if (signalling_protein == 'af-signprot') {
    segment_label = 'G Protein segments:';
  } else {
    segment_label = 'Arrestin segments:';
  };

  const headerInnerLegend = innercircleLegend.append("text")
                              .attr("x", 0)
                              .attr("y", -20)
                              .attr("class", "legend-header")
                              .text(segment_label)
                              .attr("font-family", "Arial")
                              .attr("font-size", "16px")
                              .attr("font-weight", "bold")
                              .attr("fill", "black");

  const txtInnerLen = headerInnerLegend.node().getComputedTextLength();

  let cumulativeInnerX = txtInnerLen + 20;

  Object.entries(inner_legend).forEach((type, i) => {
    const legendInnerItem = innercircleLegend.append("g")
      .attr("transform", `translate(${cumulativeInnerX}, 0)`);

    // Draw the rectangle at y=0
    legendInnerItem.append("rect")
      .attr("width", 20)
      .attr("height", 20)
      .attr("y", -35)  // <-- Set y to 0 for the rectangle
      .attr("fill", allColors[type[1]])
      .attr("stroke", "black")
      .attr("stroke-width", 1)
      .on("click", function() {
        if(d3.event) {
          d3.event.stopPropagation();
        }
        let indexValues = [];
        // Select all circles with the given segment value and gather the data-index attributes
        d3.selectAll(`circle[data-segment="${type[1]}"]`).each(function() {
          const dataIndex = d3.select(this).attr('data-index');
          if (dataIndex !== null) {
            indexValues.push(dataIndex);
          }
        });
        const sanitizedType = sanitizeClassName(type[1]);
        d3.selectAll("circle").attr("stroke-opacity", 0.2);
        d3.selectAll("path").attr("stroke-opacity", 0);
        d3.selectAll(`circle[data-segment="${type[1]}"]`).attr("stroke-opacity", 1);
        indexValues.forEach(index => {
          d3.selectAll(`path[data-inner-bead="${index}"]`).attr("stroke-opacity", 1);
        });
      });

    // Draw the text at y=15
    const textInnerElement = legendInnerItem.append("text")
      .attr("x", 25)
      .attr("y", -20)  // <-- Set y to 15 for the text
      .text(type[1]);

    // Get text width and update the cumulativeX for the next legend item
    const textInnerWidth = textInnerElement.node().getComputedTextLength();
    cumulativeInnerX += textInnerWidth + 30; // 50 is the space between the rectangle and the next item
  });

  // Creating the inner legend structure

  const headerOuterLegend = outercircleLegend.append("text")
                              .attr("x", 0)
                              .attr("y", -20)
                              .attr("class", "legend-header")
                              .text("GPCR segments:")
                              .attr("font-family", "Arial")
                              .attr("font-size", "16px")
                              .attr("font-weight", "bold")
                              .attr("fill", "black");

  const txtOuterLen = headerOuterLegend.node().getComputedTextLength();

  let cumulativeOuterX = txtOuterLen + 20;

  Object.entries(outer_legend).forEach((type, i) => {
    const legendOuterItem = outercircleLegend.append("g")
      .attr("transform", `translate(${cumulativeOuterX}, 0)`);

    // Draw the rectangle at y=0
    legendOuterItem.append("rect")
      .attr("width", 20)
      .attr("height", 20)
      .attr("y", -35)  // <-- Set y to 0 for the rectangle
      .attr("fill", allColors[type[1]])
      .attr("stroke", "black")
      .attr("stroke-width", 1)
      .on("click", function() {
        if(d3.event) {
          d3.event.stopPropagation();
        }
        let indexValues = [];
        // Select all circles with the given segment value and gather the data-index attributes
        d3.selectAll(`circle[data-segment="${type[1]}"]`).each(function() {
          const dataIndex = d3.select(this).attr('data-index');
          if (dataIndex !== null) {
            indexValues.push(dataIndex);
          }
        });
        const sanitizedType = sanitizeClassName(type[1]);
        d3.selectAll("circle").attr("stroke-opacity", 0.2);
        d3.selectAll("path").attr("stroke-opacity", 0);
        d3.selectAll(`circle[data-segment="${type[1]}"]`).attr("stroke-opacity", 1);
        indexValues.forEach(index => {
          d3.selectAll(`path[data-outer-bead="${index}"]`).attr("stroke-opacity", 1);
        });
      });

    // Draw the text at y=15
    const textOuterElement = legendOuterItem.append("text")
      .attr("x", 25)
      .attr("y", -20)  // <-- Set y to 15 for the text
      .text(type[1]);

    // Get text width and update the cumulativeX for the next legend item
    const textOuterWidth = textOuterElement.node().getComputedTextLength();
    cumulativeOuterX += textOuterWidth + 30; // 50 is the space between the rectangle and the next item
  });

  // d3.selectAll(`.hydrophobic`).attr("stroke-opacity", 0);

  function applyPresetColors(svg) {
    svg.selectAll("circle")
      .each(function(d, i) {
        const aa = d3.select(this).attr("data-aa"); // Assuming 'data-aa' stores the amino acid type
        d3.select(this).attr("fill", presetColors[aa][0]);

        // Assuming the text element immediately follows the circle
        d3.select(this.nextElementSibling).attr("fill", presetColors[aa][1]);
      });
  }

  document.getElementById("propertiesButton").addEventListener("click", function() {
    applyPresetColors(svg2); // Assuming your SVG is stored in the svg2 variable
  });

  // Function to reset colors
  function resetColors(svg) {
    svg.selectAll("circle")
      .attr("fill", "white")
      .attr("stroke-opacity", 1);
      // .attr("stroke", "black");

    svg.selectAll(`.van-der-waals`).attr("stroke-opacity", 0);

    svg.selectAll("text")
      .attr("fill", "black");

    svg.selectAll("path")
      .attr("stroke", "black")
      .attr("stroke-width", 2);

    svg.selectAll("path[segment='external']")
      .each(function(d, i) {
        const interaction = d3.select(this).attr("interaction");
        d3.select(this).attr("stroke", allColors[interaction]);
      });
  }

  // Add event listener for reset button
  document.getElementById("resetButton").addEventListener("click", function() {
    resetColors(svg2);  // Reset colors of circles
    d3.selectAll("path").attr("stroke-opacity", 1);  // Reset opacity of lines
    d3.selectAll("circle").attr("fill-opacity", 1);  // Reset opacity of lines
  });

  function applySegmentColors(svg) {
    svg.selectAll("circle[data-circle='outer']")
      .each(function(d, i) {
        const segment = d3.select(this).attr("data-segment");
        if (segment) {
          d3.select(this).attr("fill", allColors[segment]);
        }
      });
    svg.selectAll("path[segment='external']")
      .each(function(d, i) {
        const interaction = d3.select(this).attr("interaction");
        d3.select(this).attr("stroke", allColors[interaction]);
      });
  }

  document.getElementById("segmentButton").addEventListener("click", function() {
    applySegmentColors(svg2);
  });

  function applyInteractionColors(svg) {
    svg.selectAll("path")
      .filter(function() {
        return !d3.select(this).attr("segment");
      })
      .each(function(d, i) {
        const interaction = d3.select(this).attr("interaction");
        d3.select(this).attr("stroke", allColors[interaction]);
      });
  }

  document.getElementById("interactionsButton").addEventListener("click", function() {
    applyInteractionColors(svg2);
  });

}
