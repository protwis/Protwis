/*global _, d3, get_max_interface_count, xScale, bool_visible, */
// * CONSTANTS
var margin = { top: 60, right: 20, bottom: 180, left: 240 };
var w = 1200 - margin.left - margin.right,
  h = 1000 - margin.top - margin.bottom;
// change to 600 for more compact view
var non_int_col = "#fff";
// array for data in infobox
var info_data = [];
var con_seq = {};
var signprotmat = {
  // * DATA TRANSFORMATION FUNCTIONS
  data: {
    extractPdbIDs(dataset) {
      var ret = [];
      Object.keys(dataset).forEach(function (e) {
        return ret.push(e.toUpperCase());
      });
      return ret;
    },

    objectToArray(dataset) {
      return Object.keys(dataset).map(function (key) {
        return dataset[key];
      });
    },

    moveKeyToArray(dataset, pdb_ids) {
      for (var i = 0; i < pdb_ids.length; i++) {
        var pdb = pdb_ids[i];
        for (var j = 0; j < dataset[i].length; j++) {
          dataset[i][j].push(pdb);
        }
      }
      return dataset;
    },

    // https://stackoverflow.com/questions/10865025/merge-flatten-an-array-of-arrays-in-javascript/25804569#comment50580016_10865042
    flattenOnce(array) {
      return [].concat.apply([], array);
    },

    labelData(data, keys) {
      var data_labeled = data.map(function (e) {
        var obj = {};
        keys.forEach(function (key, i) {
          obj[key] = e[i];
        });
        return obj;
      });
      return data_labeled;
    },

    getInteractionTypes(dataset) {
      var int_ty = [];
      for (var i = 0; i < dataset.length; i++) {
        var int_arr = dataset[i].int_ty;
        int_ty.push(int_arr);
      }
      int_ty = signprotmat.data.flattenOnce(int_ty).filter(function (v, i, a) {
        return a.indexOf(v) === i;
      });
      var rm_index = int_ty.indexOf("undefined");
      if (rm_index > -1) {
        int_ty.splice(rm_index, 1);
      }
      return int_ty;
    },

    get_additional_receptors(data, xvals, prids) {
      var new_receptor_data = [];
      for (var index = 0; index < data.length; index++) {
        var e1 = data[index]["rec_gn"];
        var e2 = data[index]["entry_name"];
        if (xvals.includes(e1) && prids.includes(e2)) {
          new_receptor_data.push(data[index]);
        }
      }
      return new_receptor_data;
    },

    extractRecSigData(data, which_component) {
      if (which_component === "rec") {
        return _.uniqBy(data, function (t) {
          return [t.rec_gn, t.pdb_id].join();
        });
      } else if (which_component === "sig") {
        return _.uniqBy(data, function (t) {
          return [t.sig_gn, t.pdb_id].join();
        });
      } else {
        console.log("No component specified...");
      }
    },

    select_by_value(selection, value) {
      var ret_sel = [];
      for (var index = 0; index < selection.length; index++) {
        ret_sel.push(selection[index][value]);
      }
      return ret_sel;
    },

    removeUndefinedGN(dataset) {
      return _.filter(dataset, function (o) {
        return o.rec_gn !== "-";
      });
    },

    dataTransformationWrapper(dataset, pdb_sel) {
      // dataset = _.pick(dataset, pdb_sel);
      // let pdb_ids = signprotmat.data.extractPdbIDs(dataset);
      // let data_t = signprotmat.data.objectToArray(dataset);
      // data_t = signprotmat.data.moveKeyToArray(data_t, pdb_ids);
      // data_t = signprotmat.data.flattenOnce(data_t);
      // data_t = signprotmat.data.removeUndefinedGN(data_t);
      // console.log(pdb_sel)
      var data_t = _.filter(dataset, function (d) {
        return pdb_sel.includes(d.pdb_id);
      });
      var data_t_rec = signprotmat.data.extractRecSigData(data_t, "rec");
      var data_t_sig = signprotmat.data.extractRecSigData(data_t, "sig");
      var int_ty = signprotmat.data.getInteractionTypes(data_t);
      var pdb_ids = _.uniqBy(data_t, "pdb_id");
      pdb_ids = _.map(pdb_ids, function (d) {
        return d.pdb_id;
      });
      // console.log(data_t)
      // console.log(data_t_rec)
      // console.log(data_t_sig)
      // console.log(int_ty)
      // console.log(pdb_ids)
      var return_data = {
        transformed: data_t,
        receptor: data_t_rec,
        signprot: data_t_sig,
        inttypes: int_ty,
        pdbids: pdb_ids,
      };
      return return_data;
    },

    annotateNonInteractionData(meta, data) {
      data.forEach(function (element) {
        var tmp = _.find(meta, function (d) {
          return d.entry_name === element.entry_name;
        });
        element["pdb_id"] = tmp.pdb_id;
      });
      return data;
    },

    get_data_gap_non_int(pdbScale, xScale, data_non, data_receptor) {
      // Appending a gap as the background for all residues
      var data_gap = [];

      for (var pdb of pdbScale.domain()) {
        for (var gn of xScale.domain()) {
          var d = {};
          d["rec_gn"] = gn;
          d["pdb_id"] = pdb;
          d["rec_aa"] = "-";
          data_gap.push(d);
        }
      }

      data_non = _.filter(data_non, function (d) {
        return xScale(d.rec_gn);
      });

      data_non = _.filter(data_non, function (d) {
        return pdbScale(d.pdb_id);
      });

      data_gap.push(...data_non, ...data_receptor);
      // data_non.push(...data_receptor)
      // data_gap = data_non
      // console.log(data_gap)

      return data_gap;
    },

    get_receptor_sequence(data, receptor_pdb) {
      // filter all of the data for the desired receptor
      data = data.filter((data) => data.pdb_id === receptor_pdb);
      var entry_name = data.slice(-1)[0].entry_name;
      // it is assumed that the data comes from get_data_gap_non_int
      // which stacks all of the sequence information in the order of
      // first: gap for every residue
      // secon: non interacting residues the receptor has
      // last: interacting residues
      // such, reversing and uniqBy on rec_gn gives the correct sequence
      data = _.uniqBy(data.reverse(), "rec_gn");
      data = data.map((x) => ({
        rec_tm: x.rec_gn.match(/\d{1,2}/)[0],
        rec_gn: x.rec_gn.slice(-2),
        rec_aa: x.rec_aa,
      }));
      data = _.orderBy(data, ["rec_tm", "rec_gn"]);
      var seq = data.map((x) => x.rec_aa).join("");

      return entry_name.concat(">\n", seq);
    },

    combine_fasta(array_of_fasta) {
      return array_of_fasta.join("\n");
    },

    download_fasta(data, pdbScale) {
      var fasta = [];
      for (var pdb of pdbScale.domain()) {
        var entry = signprotmat.data.get_receptor_sequence(data, pdb);
        fasta.push(entry);
      }
      var fasta_export = signprotmat.data.combine_fasta(fasta);

      $("<a></a>")
        .attr("id", "downloadFile")
        .attr(
          "href",
          "data:text/fasta;charset=utf8," + encodeURIComponent(fasta_export)
        )
        .attr("download", "export.fasta")
        .appendTo("body");

      $("#downloadFile").ready(function () {
        $("#downloadFile").get(0).click();
        $("#downloadFile").remove();
      });
    },
  },

  // * D3 DRAW FUNCTIONS
  d3: {
    // * SETTING UP SVG FOR OUTPUT
    setup(div, loc) {
      if (loc === "seqsig") {
        h = 1200;
        var mt = 0;
      } else if (loc === "conseq") {
        h = 0;
        var mt = 30;
      } else {
        loc = "interaction";
        h = 1000 - margin.top - margin.bottom;
        var mt = 100;
      }
      var svg = d3
        .select("body")
        .select("div#content")
        .select(div)
        .append("svg")
        .attr("preserveAspectRatio", "xMinYMin meet")
        .attr(
          "viewBox",
          "0 0 " +
            (w + margin.left + margin.right) +
            " " +
            (h + 200 + margin.top + margin.bottom)
        )
        // .classed("svg-content", true) //class to make it responsive
        .attr(
          "class",
          typeof loc !== "undefined" ? "svg-content " + loc : "svg-content"
        )
        .append("g")
        .attr("transform", "translate(" + margin.left + "," + mt + ")");
      return svg;
    },

    // * SETTING THE X/Y SCALE
    xScale(data, receptor) {
      var domain = d3
        .map(data, function (d) {
          return d.rec_gn;
        })
        .keys()
        .sort(function (a, b) {
          var a_patt = /(\d*)x/g;
          var b_patt = /(\d*)x/g;
          var a_match = a_patt.exec(a);
          var b_match = b_patt.exec(b);
          var a_obj = _.findIndex(receptor, function (d) {
            if (a_match){
              return d === a_match[1];
            } else if (a === '-'){
              return d === a;
            }
          });
          var b_obj = _.findIndex(receptor, function (d) {
            if (b_match){
              return d === b_match[1];
            } else if (b === '-'){
              return d === b;
            }
          });
          // console.log(a_obj);
          return d3.ascending(a_obj, b_obj);
        });
      var xScale = d3
        .scaleBand()
        .domain(domain)
        .range([0, w])
        // .round(true)
        .padding(1);
      return xScale;
    },

    yScale(data, gprot) {
      var domain = d3
        .map(data, function (d) {
          return d.sig_gn;
        })
        .keys()
        // .sort(d3.descending);
        .sort(function (a, b) {
          // TEMP FIX for Arrestins - when N. and C. tag is added to segment name, this should be removed
          if (a.startsWith('N.') || a.startsWith('C.')) {
            var a_patt = /\.(\S*)\./g;
            var b_patt = /\.(\S*)\./g;
          }
          else {
            var a_patt = /(\S*)\./g;
            var b_patt = /(\S*)\./g;
          }

          var a_match = a_patt.exec(a);
          var b_match = b_patt.exec(b);
          var a_obj = _.find(gprot, function (d) {
            return d.slug === a_match[1];
          });
          var b_obj = _.find(gprot, function (d) {
            return d.slug === b_match[1];
          });
          // console.log(a_obj.id);
          // console.log(a_obj.id + (+a.slice(-2)/100));
          return d3.descending(
            a_obj.id + +a.slice(-2) / 100,
            b_obj.id + +b.slice(-2) / 100
          );
          // return d3.descending(a_obj.id, b_obj.id);
        });
      var yScale = d3.scaleBand().domain(domain).range([h, 0]).padding(1);
      return yScale;
    },

    // * SETTING THE PDB/SIG-PROT SCALE
    pdbScale(data, meta) {
      var pdbScale = d3
        .scaleBand()
        .domain(
          d3
            .map(data, function (d) {
              return d.pdb_id;
            })
            .keys()
            .sort(function (a, b) {
              var a_obj = _.find(meta, function (d) {
                return d.pdb_id === a;
              });
              var b_obj = _.find(meta, function (d) {
                return d.pdb_id === b;
              });
              return d3.descending(a_obj.entry_name, b_obj.entry_name);
            })
        )
        .range([300, 0])
        .padding(1);
      return pdbScale;
    },

    sigScale(data, meta) {
      // calculating the number of columns in the scale
      var pdbs = [];
      var range;
      for (var i=0; i < data.length; i++){
        if (!pdbs.includes(data[i]["pdb_id"])){
          pdbs.push(data[i]["pdb_id"]);
        }
      }
      if (pdbs.length <= 300){
        range = 120;
      } else {
        range = 12*pdbs.length;
      }
      var sigScale = d3
        .scaleBand()
        .domain(
          d3
            .map(data, function (d) {
              return d.pdb_id;
            })
            .keys()
            .sort(function (a, b) {
              var a_obj = _.find(meta, function (d) {
                return d.pdb_id === a;
              });
              var b_obj = _.find(meta, function (d) {
                return d.pdb_id === b;
              });
              return d3.descending(a_obj.gprot, b_obj.gprot);
            })
        )
        .range([range, 0])
        .padding(1);
      return sigScale;
    },

    // * SETTING THE COLOR SCALE
    colScale(f) {
      var scale = {
        "van-der-waals": "#d9d9d9",
        // "edge-to-face": "#969696",
        // "water-mediated": "#7DB144",
        hydrophobic: "#A6E15F", //"#93d050",
        // "polar-sidechain-sidechain": "#EAA91F",
        // "polar-sidechain-backbone": "#C38E1A",
        // "polar-backbone-sidechain": "#C3A563",
        // "h-bond donor-acceptor": "#7030a0",
        // "h-bond acceptor-donor": "#B24DFF",
        // "cation-pi": "#0070c0",
        // "pi-cation": "#005693",
        ionic: "#005693",
        aromatic: "#689235", //"#7DB144",
        polar: "#7030a0",
      };
      var colScale = d3
        .scaleOrdinal()
        .domain(Object.keys(scale))
        .range(Object.values(scale));
      return colScale;
    },

    // * seqsig
    // * SETTING THE FEATURE SCALE
    fScale(data) {
      var f = _.map(data, function (d) {
        var length_text = d.length != "" ? " (" + d.length + ")" : "";
        return d.feature + length_text;
      });
      var fScale = d3
        .scaleBand()
        .domain(f)
        .range([0, h])
        // .round(true)
        .padding(1);

      return fScale;
    },

    resScaleColor(f) {
      var scale = {
        A: {bg_color: "#B3B300", font_color: "#000000"},
        C: {bg_color: "#7B7D3A", font_color: "#000000"},
        D: {bg_color: "#A50707", font_color: "#FDFF7B"},
        E: {bg_color: "#A50707", font_color: "#FDFF7B"},
        F: {bg_color: "#11B90A", font_color: "#000000"},
        G: {bg_color: "#FF00C8", font_color: "#000000"},
        H: {bg_color: "#0075A5", font_color: "#000000"},
        I: {bg_color: "#B3B300", font_color: "#000000"},
        K: {bg_color: "#0B3E9E", font_color: "#FDFF7B"},
        L: {bg_color: "#B3B300", font_color: "#000000"},
        M: {bg_color: "#B3B300", font_color: "#000000"},
        N: {bg_color: "#880896", font_color: "#FDFF7B"},
        P: {bg_color: "#990066", font_color: "#FDFF7B"},
        Q: {bg_color: "#880896", font_color: "#FDFF7B"},
        R: {bg_color: "#0B3E9E", font_color: "#FDFF7B"},
        S: {bg_color: "#880896", font_color: "#FDFF7B"},
        T: {bg_color: "#880896", font_color: "#FDFF7B"},
        V: {bg_color: "#B3B300", font_color: "#000000"},
        W: {bg_color: "#0B7C00", font_color: "#000000"},
        Y: {bg_color: "#11B90A", font_color: "#000000"},
        "-": {bg_color: "#FFFFFF", font_color: "#000000"},
        _: {bg_color: "#EDEDED", font_color: "#000000"},
        "+": {bg_color: "#FFFFFF", font_color: "#000000"},
      };

      return scale[f];
    },

    print_resScaleColor_legend() {
      var colScale = signprotmat.d3.colScale();
      var scale = {
        A: {bg_color: "#B3B300", font_color: "#000000"},
        C: {bg_color: "#7B7D3A", font_color: "#000000"},
        D: {bg_color: "#A50707", font_color: "#FDFF7B"},
        E: {bg_color: "#A50707", font_color: "#FDFF7B"},
        F: {bg_color: "#11B90A", font_color: "#000000"},
        G: {bg_color: "#FF00C8", font_color: "#000000"},
        H: {bg_color: "#0075A5", font_color: "#000000"},
        I: {bg_color: "#B3B300", font_color: "#000000"},
        K: {bg_color: "#0B3E9E", font_color: "#FDFF7B"},
        L: {bg_color: "#B3B300", font_color: "#000000"},
        M: {bg_color: "#B3B300", font_color: "#000000"},
        N: {bg_color: "#880896", font_color: "#FDFF7B"},
        P: {bg_color: "#990066", font_color: "#FDFF7B"},
        Q: {bg_color: "#880896", font_color: "#FDFF7B"},
        R: {bg_color: "#0B3E9E", font_color: "#FDFF7B"},
        S: {bg_color: "#880896", font_color: "#FDFF7B"},
        T: {bg_color: "#880896", font_color: "#FDFF7B"},
        V: {bg_color: "#B3B300", font_color: "#000000"},
        W: {bg_color: "#0B7C00", font_color: "#000000"},
        Y: {bg_color: "#11B90A", font_color: "#000000"},
        "-": {bg_color: "#FFFFFF", font_color: "#000000"},
        _: {bg_color: "#EDEDED", font_color: "#000000"},
        "+": {bg_color: "#FFFFFF", font_color: "#000000"},
      };

      function responsivefy(svg) {
        // get container + svg aspect ratio
        var container = d3.select(svg.node().parentNode),
          width = parseInt(svg.style("width")),
          height = parseInt(svg.style("height")),
          aspect = width / height;

        // add viewBox and preserveAspectRatio properties,
        // and call resize so that svg resizes on inital page load
        svg
          .attr("viewBox", "0 0 " + width + " " + height)
          .attr("perserveAspectRatio", "xMinYMid")
          .call(resize);

        // to register multiple listeners for same event type,
        // you need to add namespace, i.e., 'click.foo'
        // necessary if you call invoke this function for multiple svgs
        // api docs: https://github.com/mbostock/d3/wiki/Selections#on
        d3.select(window).on("resize." + container.attr("id"), resize);

        // get width of container and resize svg to fit it
        function resize() {
          var targetWidth = parseInt(container.style("width"));
          svg.attr("width", targetWidth);
          svg.attr("height", Math.round(targetWidth / aspect));
        }
      }

      var ordinal = d3
        .scaleOrdinal()
        .domain(Object.keys(scale))
        .range(Object.values(scale).map((x) => x.bg_color));

      var svg = d3
        .select("#legend-space")
        .append("svg")
        .attr("width", 554)
        .attr("height", 110)
        .call(responsivefy);

      svg
        .append("g")
        .attr("class", "legendOrdinal")
        .attr("transform", "translate(10,20)");

      var legendOrdinal = d3
        .legendColor()
        .orient("horizontal")
        .labelAlign("center")
        .shapePadding(2)
        .scale(ordinal);

      svg
        .select(".legendOrdinal")
        .call(legendOrdinal)
        .selectAll("rect")
        .attr("rx", 3)
        .attr("ry", 3);
      svg.select(".legendOrdinal").selectAll("text").attr("class", "legend");

      // * ADDING Interaction Type LEGEND
      let size = 2;
      let window_starts = _.range(0, colScale.domain().length + 1, size);

      let i = 0;
      for (let windo of window_starts) {
        let start = windo;
        let stop = windo + size;
        let element_ids = _.range(start, stop);
        let filter_elements = _.pullAt(colScale.domain(), element_ids);

        svg
          .append("g")
          .attr("class", "legendOrdinal" + i)
          .attr(
            "transform",
            "translate(" +
              // (xScale.step() / 2 + i * 10 * xScale.step()) + ","
              (10 + i * 160) +
              "," +
              65 +
              ")"
          );

        let legendOrdinal = d3
          .legendColor()
          .cellFilter(function (d) {
            return filter_elements.includes(d.label);
          })
          .orient("vertical")
          .labelAlign("start")
          .shapePadding(2)
          .scale(colScale);
        svg
          .select(".legendOrdinal" + i)
          .call(legendOrdinal)
          .selectAll("rect")
          .attr("rx", 3)
          .attr("ry", 3);
        svg
          .select(".legendOrdinal" + i)
          .selectAll("text")
          .attr("class", "legend");

        i += 1;
      }
    },

    fScaleColor(f) {
      if (f === "αH") {
        f = "aH";
      }
      var scale = {
        HY: { bg_color: "#93d050" },
        HA: { bg_color: "#ffff00" },
        M: { bg_color: "#ffff00" },
        A: { bg_color: "#ffff00" },
        I: { bg_color: "#ffff00" },
        L: { bg_color: "#ffff00" },
        V: { bg_color: "#ffff00" },
        HR: { bg_color: "#07b050" },
        W: { bg_color: "#07b050" },
        Y: { bg_color: "#07b050" },
        F: { bg_color: "#07b050" },
        Hb: { bg_color: "#7030a0", font_color: "#ffffff" },
        N: { bg_color: "#7030a0", font_color: "#ffffff" },
        Q: { bg_color: "#7030a0", font_color: "#ffffff" },
        S: { bg_color: "#7030a0", font_color: "#ffffff" },
        T: { bg_color: "#7030a0", font_color: "#ffffff" },
        Hu: { bg_color: "#7030a0", font_color: "#ffffff" },
        Ha: { bg_color: "#7030a0", font_color: "#ff0000" },
        Hd: { bg_color: "#7030a0", font_color: "#02b0f0" },
        "+-": { bg_color: "#0070c0", font_color: "#ff0000" },
        "+": { bg_color: "#0070c0", font_color: "#000000" },
        H: { bg_color: "#0070c0", font_color: "#000000" },
        K: { bg_color: "#0070c0", font_color: "#000000" },
        R: { bg_color: "#0070c0", font_color: "#000000" },
        "-": { bg_color: "#ff0000" },
        D: { bg_color: "#ff0000" },
        E: { bg_color: "#ff0000" },
        Sm: { bg_color: "#ffffff" },
        aH: { bg_color: "#d9d9d9" },
        G: { bg_color: "#ff02ff" },
        P: { bg_color: "#d603ff", font_color: "#ffffff" },
        C: { bg_color: "#bf8f00" },
      };
      return scale[f];
    },

    colorBySwitch(which, colScale) {
      if (which === "res") {
        var other = "int";
        var res = d3.select("g#recAA").selectAll("g");
        res.selectAll("rect").style("fill", function (d) {
          var col = signprotmat.d3.resScaleColor(d.rec_aa);
          if (typeof col !== "undefined") {
            return col.bg_color;
          } else {
            return null;
          }
        });
        res.selectAll("text").style("fill", function (d) {
          var col = signprotmat.d3.resScaleColor(d.rec_aa);
          if (typeof col !== "undefined") {
            if (typeof col.font_color !== "undefined") {
              return col.font_color;
            } else {
              return "#000000";
            }
          } else {
            return "#000000";
          }
        });
        res = d3.select("g#sigAA").selectAll("g");
        res.selectAll("rect").style("fill", function (d) {
          var col = signprotmat.d3.resScaleColor(d.sig_aa);
          if (typeof col !== "undefined") {
            return col.bg_color;
          } else {
            return null;
          }
        });
        res.selectAll("text").style("fill", function (d) {
          var col = signprotmat.d3.resScaleColor(d.sig_aa);
          if (typeof col !== "undefined") {
            if (typeof col.font_color !== "undefined") {
              return col.font_color;
            } else {
              return "#000000";
            }
          } else {
            return "#000000";
          }
        });
      } else if (which === "int") {
        var other = "res";
        var res = d3.select("g#recAA").selectAll("g");
        res.selectAll("rect").style("fill", function (d) {
          if (typeof d.int_ty !== "undefined") {
            return colScale(d.int_ty[0]);
          } else {
            return non_int_col;
          }
        });
        res.selectAll("text").style("fill", null);
        res = d3.select("g#sigAA").selectAll("g");
        res.selectAll("rect").style("fill", function (d) {
          if (typeof d.int_ty !== "undefined") {
            return colScale(d.int_ty[0]);
          } else {
            return "#ccc";
          }
        });
        res.selectAll("text").style("fill", null);
      }
      document.querySelector("#" + which + "but").classList.add("active");
      document.querySelector("#" + other + "but").classList.remove("active");
    },

    cScale(data) {
      // const values = d3
      //   .map(data, (d: any) => d.cons)
      //   .keys()
      //   .map(Number)
      //   .sort(d3.ascending)
      // const min = values.pop()
      // const max = values[0]
      // conservation is calculated to be between -1 and 10 by python
      var cScale = d3.scaleSequential(d3.interpolateGreys).domain([0, 100]);
      return cScale;
    },

    // * DEFINING AXIS FOR X/Y AND GRID
    xAxis(xScale) {
      var xAxis = d3.axisBottom(xScale).tickSize(0).tickPadding(8);
      return xAxis;
    },

    yAxis(yScale) {
      var yAxis = d3.axisLeft(yScale).tickSize(0).tickPadding(8);
      return yAxis;
    },

    xAxisGrid(xScale, yScale) {
      var xAxisGrid = d3
        .axisTop(xScale)
        .tickSize(h - yScale.step())
        .tickFormat(function (d) {
          return "";
        });
      return xAxisGrid;
    },

    yAxisGrid(xScale, yScale) {
      var yAxisGrid = d3
        .axisRight(yScale)
        .tickSize(w - xScale.step())
        .tickFormat(function (d) {
          return "";
        });
      return yAxisGrid;
    },

    // * ADD TOOLTIP FUNCTIONALITY
    tooltip(svg) {
      var tip = d3
        .tip()
        .attr("class", "d3-tip")
        .html(function (d) {
          var pair_string = "";
          d.pairs.forEach(function (e) {
            pair_string +=
              e.pdb_id +
              ": " +
              e.rec_aa +
              " vs. " +
              e.sig_aa +
              "(" +
              e.int_ty +
              ")" +
              "<br>";
          });
          return (
            "Receptor: " +
            d.rec_gn +
            "<br>" +
            "Signaling Protein: " +
            d.sig_gn +
            "<br>" +
            "PDBs: " +
            "<br>" +
            pair_string
          );
        });
      svg.call(tip);
      return tip;
    },

    // * RENDER DATA
    renderData(
      svg,
      data,
      data_non,
      interactions_metadata,
      xScale,
      yScale,
      xAxis,
      yAxis,
      xAxisGrid,
      yAxisGrid,
      colScale,
      pdbScale,
      sigScale,
      tip
    ) {
      var shift_left = 7 / 8;
      var shift_top = 1 / 8;
      var scale_size = shift_left - shift_top;
      var offset = 1;
      var each_res;
      var helper = {};
      var result = data.transformed.reduce(function (r, o) {
        var key = o.rec_gn + "-" + o.sig_gn;
        if (!helper[key]) {
          var tmp = {
            rec_gn: o.rec_gn,
            sig_gn: o.sig_gn,
            pairs: [
              {
                pdb_id: o.pdb_id,
                rec_aa: o.rec_aa,
                sig_aa: o.sig_aa,
                int_ty: o.int_ty,
              },
            ],
          };
          helper[key] = tmp;
          r.push(helper[key]);
        } else {
          helper[key].pairs.push({
            pdb_id: o.pdb_id,
            rec_aa: o.rec_aa,
            sig_aa: o.sig_aa,
            int_ty: o.int_ty,
          });
        }
        return r;
      }, []);

      var bwScale = d3
        .scaleSequential(d3.interpolateGreys)
        .domain([0, pdbScale.domain().length]);
      svg.append("g").attr("id", "interact");

      var each_rect = svg
        .select("g#interact")
        .selectAll("rects")
        .data(result)
        .enter()
        .append("g")
        .attr("class", function (d) {
          return "p" + d.pairs.length;
        })
        .on("mouseover", function (d) {
          tip.show(d);
        })
        .on("mouseout", function (d) {
          tip.hide();
        });
      // .on("click", function (d) {
      // var index;
      // // let rect_x = d3.event.target.getAttribute('x')
      // // let rect_y = d3.event.target.getAttribute('y')
      // // console.log(rect_x, rect_y)
      // // https://stackoverflow.com/a/20251369/8160230
      // // select the rect under cursor
      // var curr = d3.select(this);
      // // Determine if current rect was clicked before
      // var active = d.active ? false : true;
      // // Update whether or not the elements are active
      // d.active = active;
      // // set style in regards to active
      // if (d.active) {
      //     curr.style("stroke", "yellow").style("stroke-width", 2);
      //     info_data.push(d);
      // }
      // else {
      //     curr.style("stroke", "none").style("stroke-width", 2);
      //     index = info_data.indexOf(d);
      //     info_data.splice(index, 1);
      // }
      // signprotmat.d3.infoBoxUpdate();
      // signprotmat.d3.colorRecResidues(d);
      // signprotmat.d3.colorSigResidues(d);
      // });

      each_rect
        .append("rect")
        .attr("x", function (d) {
          // return xScale(d.rec_gn) - shift_left * xScale.step() + offset;
          return xScale(d.rec_gn) - xScale.step();
        })
        .attr("y", function (d) {
          // return yScale(d.sig_gn) + shift_top * yScale.step() + offset;
          return yScale(d.sig_gn);
        })
        .attr("rx", function () {
          if (data.transformed.length < 15) {
            return 5;
          } else {
            return 3;
          }
        })
        .attr("ry", function () {
          if (data.transformed.length < 15) {
            return 5;
          } else {
            return 3;
          }
        })
        // .attr("width", xScale.step() * scale_size)
        // .attr("height", yScale.step() * scale_size)
        .attr("width", xScale.step())
        .attr("height", yScale.step())
        .attr("fill", function (d) {
          if (pdbScale.domain().length > 1) {
            return bwScale(d.pairs.length);
          } else {
            return colScale(d.pairs[0].int_ty[0]);
          }
        });

      each_rect
        .append("text")
        .attr("text-anchor", "middle")
        .attr("class", "int_count")
        .attr("x", function (d) {
          // return xScale(d.rec_gn) - shift_left * xScale.step() + offset;
          return xScale(d.rec_gn) - xScale.step();
        })
        .attr("y", function (d) {
          // return yScale(d.sig_gn) + shift_top * yScale.step() + offset;
          return yScale(d.sig_gn);
        })
        .attr("dx", xScale.step() / 2)
        .attr("dy", yScale.step() / 2)
        .style("fill", function (d) {
          const num_pairs = d.pairs.length;
          const max_count = get_max_interface_count();
          const ratio = (num_pairs / max_count) * 100;
          if (ratio >= 50) {
            return "#eaeaea";
          } else if (ratio < 50) {
            return "#000000";
          }
        })
        // .text(function(d){ return d.pairs.length })
        .text(function (d) {
          const num_pairs = d.pairs.length;
          const max_count = get_max_interface_count();
          const ratio = (num_pairs / max_count) * 100;
          if (ratio > 100){
            return "";
          } else {
            return _.round(ratio, 0);
          }
        });

      // * DRAWING AXES
      svg
        .append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(" + -xScale.step() / 2 + "," + h + ")")
        .call(xAxis)
        .selectAll("text")
        .attr("text-anchor", "end")
        .attr("font-size", "12px")
        .attr("dx", "-5px")
        .attr("dy", "-5px")
        .attr("transform", "rotate(-90)");
      svg
        .append("g")
        .attr("class", "y axis")
        .attr("transform", "translate(" + 0 + "," + yScale.step() / 2 + ")")
        .call(yAxis)
        .selectAll("text")
        .attr("font-size", "12px");
      // * DRAWING GRIDLINES
      svg
        .append("g")
        .attr("class", "x grid")
        .attr("transform", "translate(" + 0 + "," + h + ")")
        .call(xAxisGrid);
      svg
        .append("g")
        .attr("class", "y grid")
        .attr("transform", "translate(" + 0 + "," + yScale.step() + ")")
        .call(yAxisGrid);
      // * ADDITIONAL FIGURE LINES
      // top x line
      svg
        .append("line")
        .style("stroke", "black")
        .attr("x1", 0)
        .attr("y1", yScale.step())
        .attr("x2", w - xScale.step())
        .attr("y2", yScale.step());
      // left y line
      svg
        .append("line")
        .style("stroke", "black")
        .attr("x1", w - xScale.step())
        .attr("y1", yScale.step())
        .attr("x2", w - xScale.step())
        .attr("y2", h);
      // * ADD AXIS LABELS
      svg
        .append("text")
        .attr("class", "x axis_label")
        .attr("text-anchor", "end")
        .attr("x", 0 - 7)
        .attr("y", h + 25)
        .text("Res. Pos.");
      // svg
      //     .append("text")
      //     .attr("class", "y axis_label")
      //     .attr("text-anchor", "end")
      //     .attr("x", 0 - 7)
      //     .attr("y", 0.8 * yScale.step())
      //     .text("Res. Pos.");
      // * ADD INFOBOX ELEMENT
      svg
        .append("g")
        .attr("id", "infobox")
        .attr(
          "transform",
          "translate(-15," + (data.inttypes.length + 2) * 20 + ")"
        );

      // * APPENDING COL TICK ANNOTATION FOR RECEPTOR GNs
      let greek_match = {
        "&alpha;": "\u03B1",
        "&beta;": "\u03B2",
        "&gamma;": "\u03B3",
        "&delta;": "\u03B4",
        "&kappa;": "\u03BA",
        "&mu;": "\u03BC",
      };
      svg
        .append("g")
        .attr("id", "recPDB")
        .attr("transform", "translate(" + 0 + "," + h * 0.97 + ")")
        .selectAll("text")
        .data(data.pdbids)
        .enter()
        .append("text")
        .attr("class", "y seq_label")
        .attr("x", -10)
        .attr("y", function (d) {
          return pdbScale(d) - pdbScale.step() / 2;
        })
        .attr("text-anchor", "end")
        .attr("dy", 75)
        .text(function (d) {
          var i_obj = _.find(interactions_metadata, function (e) {
            return e.pdb_id === d;
          });
          var text = i_obj.name;
          // replace greek letters
          for (var code in greek_match) {
            if (Object.prototype.hasOwnProperty.call(greek_match, code)) {
              text = text.replace(code, greek_match[code]);
            }
          }
          return text.replace(/<[^>]*>/g, "") + " (" + d.toUpperCase() + ")";
          // return d;
        });
      // * APPENDING ROW TICK ANNOTATION FOR SIGPROT GNs
      svg
        .append("g")
        .attr("id", "sigPDB")
        .attr(
          "transform",
          "translate(" +
            (0 - margin.left * 0.75 - sigScale.range()[0] / 2 - 5) +
            "," +
            yScale.step() +
            ")rotate(-90)"
        )
        .selectAll("text")
        .data(data.pdbids)
        .enter()
        .append("text")
        .attr("class", "x seq_label")
        .attr("x", function (d, i) {
          return 10;
        })
        .attr("y", function (d, i) {
          // return sigScale(d) - sigScale.step() / 2;
          return sigScale(d);
          // return sigScale.step() * (i + 1);
        })
        .attr("text-anchor", "begin")
        .attr("dy", 68)
        .text(function (d) {
          var i_obj = _.find(interactions_metadata, function (e) {
            return e.pdb_id === d;
          });
          // let text = i_obj.gprot.replace('Engineered', 'Eng.')
          var text = i_obj.gprot.replace("Engineered", "E.");
          // text = text.replace('protein', 'prot.')
          text = text.replace("protein", "p.");
          return text.replace(/<[^>]*>/g, "") + " (" + d.toUpperCase() + ")";
        });
      // * APPENDING AMINOACID SEQUENCE [RECEPTOR]
      var recTip = d3
        .tip()
        .attr("class", "d3-tip")
        .html(function (d) {
          return (
            "Signal Prot. AA: " +
            d.sig_aa +
            "<br>" +
            "Interaction type: " +
            d.int_ty
          );
        });
      svg
        .append("g")
        .attr("id", "recAA")
        .attr(
          "transform",
          "translate(" + -xScale.step() / 2 + "," + h * 0.97 + ")"
        )
        .append("rect")
        .attr("class", "border-bg")
        .style("fill", "#ffffff")
        .attr("x", xScale.step() / 2)
        .attr("y", 75)
        .attr("width", xScale.range()[1] - xScale.step())
        .attr("height", pdbScale.range()[0] - pdbScale.step());

      var data_gap = signprotmat.data.get_data_gap_non_int(
        pdbScale,
        xScale,
        data_non,
        data.receptor
      );

      //Appending residues that are non interacting and residues that do interact on top
      each_res = svg
        .select("g#recAA")
        .selectAll("text")
        .data(data_gap)
        .enter()
        .append("g")
        .attr("class", function (d) {
          return "R_" + _.replace(d.rec_gn, ".", "p") + "_P_" + d.pdb_id;
        })
        .call(recTip)
        .on("mouseover", function (d) {
          if (d.int_ty != null) {
            recTip.show(d);
          }
        })
        .on("mouseout", function (d) {
          recTip.hide();
        });

      each_res
        .append("rect")
        .attr("class", "res_rect")
        .style("fill", function (d) {
          var col = signprotmat.d3.resScaleColor(d.rec_aa); // Assuming d.rec_aa is the residue data for receptor
          if (typeof col !== "undefined") {
            return col.bg_color;
          } else {
            return null; // Or a default color if undefined
          }
        })
        .attr("x", function (d) {
          return xScale(d.rec_gn) - xScale.step() / 2;
        })
        .attr("y", function (d) {
          return 75 + pdbScale(d.pdb_id) - pdbScale.step();
        })
        .attr("width", xScale.step())
        .attr("height", pdbScale.step());
      each_res
        .append("text")
        .attr("class", "res_label")
        .attr("x", function (d) {
          return xScale(d.rec_gn);
        })
        .attr("y", function (d) {
          return pdbScale(d.pdb_id) - pdbScale.step() / 2;
        })
        .attr("text-anchor", "middle")
        .attr("dy", 75)
        .style("fill", function (d) {
          return typeof d.int_ty !== "undefined"
            ? non_int_col
            : "black";
        })
        .text(function (d) {
          return d.rec_aa;
        });
      // d3.select("g#recAA")
      //     .append("rect")
      //     .attr("class", "border")
      //     .style("stroke", "black")
      //     .style("fill", "none")
      //     .attr("x", xScale.step() / 2)
      //     .attr("y", 75)
      //     .attr("width", xScale.range()[1] - xScale.step())
      //     .attr("height", pdbScale.range()[0] - pdbScale.step());

      // * APPENDING AMINOACID SEQUENCE [SIGPROT]
      var sigTip = d3
        .tip()
        .attr("class", "d3-tip")
        .html(function (d) {
          return (
            "Receptor AA: " +
            d.rec_aa +
            "<br>" +
            "Interaction type: " +
            d.int_ty
          );
        });
      svg
        .append("g")
        .attr("id", "sigAA")
        .attr(
          "transform",
          "translate(" +
            // (w + (1 / 3) * margin.right) +
            (0 - margin.left * 0.75) +
            "," +
            yScale.step() / 2 +
            ")"
        )
        .append("rect")
        .style("fill", "#ffffff")
        .attr("x", 0 + sigScale.step() / 2)
        .attr("y", yScale.step() / 2)
        .attr("width", sigScale.range()[0] - sigScale.step())
        .attr("height", yScale.range()[0] - yScale.step());
      each_res = svg
        .select("g#sigAA")
        .selectAll("text")
        .data(data.signprot)
        .enter()
        .append("g")
        .attr("class", function (d) {
          return "S_" + d.sig_gn.replace(/\./gi, "p") + "_P_" + d.pdb_id;
        })
        .call(sigTip)
        .on("mouseover", function (d) {
          if (d.int_ty !== "undefined") {
            sigTip.show(d);
          }
        })
        .on("mouseout", function (d) {
          sigTip.hide();
        });
      each_res
        .append("rect")
        .attr("class", "res_rect_vertical")
        .style("fill", function (d) {
          var col = signprotmat.d3.resScaleColor(d.sig_aa); // Assuming d.rec_aa is the residue data for receptor
          if (typeof col !== "undefined") {
            return col.bg_color;
          } else {
            return null; // Or a default color if undefined
          }
        })
        .attr("x", function (d) {
          return sigScale(d.pdb_id) - sigScale.step() / 2;
        })
        .attr("y", function (d) {
          return yScale(d.sig_gn) - yScale.step() / 2;
        })
        .attr("width", sigScale.step())
        .attr("height", yScale.step());
      each_res
        .append("text")
        .attr("class", "res_label")
        .style("fill", "white")
        .attr("x", function (d) {
          return sigScale(d.pdb_id);
        })
        .attr("y", function (d) {
          return yScale(d.sig_gn);
        })
        .attr("text-anchor", "middle")
        .text(function (d) {
          return d.sig_aa;
        });
      // d3.select("g#sigAA")
      //     .append("rect")
      //     .attr("class", "border")
      //     .style("stroke", "black")
      //     .style("fill", "none")
      //     .attr("x", 0 + sigScale.step() / 2)
      //     .attr("y", yScale.step() / 2)
      //     .attr("width", sigScale.range()[0] - sigScale.step())
      //     .attr("height", yScale.range()[0] - yScale.step());
      return svg;
    },

    addReceptor(new_data, data, svg) {
      data = _.union(data.transformed, new_data);
      data = signprotmat.data.extractRecSigData(data, "rec");
      var pdb_ids = [];
      _.forEach(_.uniqBy(data, "pdb_id"), function (value) {
        pdb_ids.push(value["pdb_id"]);
      });
      var pdbScale = signprotmat.d3.pdbScale(data);
      var xScale = signprotmat.d3.xScale(data);
      var selection = svg
        .select("g#recAA")
        .selectAll("text.res_label")
        .data(data);
      var selection_rect = svg
        .select("g#recAA")
        .selectAll("rect.res_rect")
        .data(data);
      var selection_enter = selection.enter().append("g");
      selection_enter
        .append("rect")
        .attr("class", "res_rect")
        .style("fill", "slategrey")
        .attr("x", function (d) {
          return xScale(d.rec_gn) - xScale.step() / 2;
        })
        .attr("y", function (d) {
          return 75 + pdbScale(d.pdb_id) - pdbScale.step();
        })
        .attr("width", xScale.step())
        .attr("height", pdbScale.step())
        .merge(selection_rect)
        .transition()
        .duration(500)
        .attr("x", function (d) {
          return xScale(d.rec_gn) - xScale.step() / 2;
        })
        .attr("y", function (d) {
          return 75 + pdbScale(d.pdb_id) - pdbScale.step();
        })
        .attr("width", xScale.step())
        .attr("height", pdbScale.step());
      selection_enter
        .append("text")
        .attr("class", "res_label")
        .style("fill", "white")
        .attr("x", function (d) {
          return xScale(d.rec_gn);
        })
        .attr("y", function (d) {
          return pdbScale(d.pdb_id) - pdbScale.step() / 2;
        })
        .attr("text-anchor", "middle")
        .attr("dy", 75)
        .text(function (d) {
          return d.rec_aa;
        })
        .merge(selection)
        .transition()
        .duration(500)
        .attr("x", function (d) {
          return xScale(d.rec_gn);
        })
        .attr("y", function (d) {
          return pdbScale(d.pdb_id) - pdbScale.step() / 2;
        });
      selection.exit().transition().duration(500).remove();
      selection_rect.exit().transition().duration(500).remove();
      selection = svg.select("g#recPDB").selectAll("text").data(pdb_ids);
      selection_enter = selection.enter();
      selection_enter
        .append("text")
        .attr("class", "y seq_label")
        .attr("x", -10)
        .attr("y", function (d) {
          return pdbScale(d) - pdbScale.step() / 2;
        })
        .attr("text-anchor", "end")
        .attr("dy", 75)
        .text(function (d) {
          return d;
        })
        .merge(selection)
        .transition()
        .duration(500)
        .attr("x", -10)
        .attr("y", function (d) {
          return pdbScale(d) - pdbScale.step() / 2;
        })
        .attr("dy", 75);
      selection.exit().transition().duration(500).remove();
      d3.select("g#recAA")
        .selectAll("rect.border, rect.border-bg")
        .transition()
        .duration(500)
        .attr("x", xScale.step() / 2)
        .attr("y", 75)
        .attr("width", xScale.range()[1] - xScale.step())
        .attr("height", pdbScale.range()[0] - pdbScale.step());
    },

    colorRecResidues(d) {
      var rec_gn = _.replace(d.rec_gn, ".", "p");
      var pdb_list = d.pairs.map(function (x) {
        return x["pdb_id"];
      });
      // select the rect in the g that corresponds to this rec_gn and pdb_id
      pdb_list.forEach(function (pdb) {
        d3.select("g." + "R_" + rec_gn + "_P_" + pdb)
          .select("rect")
          .classed("activeRes", d.active ? true : false);
      });
    },

    colorSigResidues(d) {
      var sig_gn = d.sig_gn.replace(/\./gi, "p");
      var pdb_list = d.pairs.map(function (x) {
        return x["pdb_id"];
      });
      // select the rect in the g that corresponds to this rec_gn and pdb_id
      pdb_list.forEach(function (pdb) {
        d3.select("g." + "S_" + sig_gn + "_P_" + pdb)
          .select("rect")
          .classed("activeRes", d.active ? true : false);
      });
    },

    infoBoxUpdate() {
      // create selection and bind data
      var info_box = d3.select("g#infobox").selectAll("text").data(info_data);
      // update existing nodes
      info_box
        .attr("y", function (d, i) {
          return i * 15;
        })
        .attr("text-anchor", "end")
        .attr("class", "legend");
      // create nodes for new data
      info_box
        .enter()
        .append("text")
        .attr("y", function (d, i) {
          return i * 15;
        })
        .attr("text-anchor", "end")
        .attr("class", "legend")
        .text(function (d) {
          return d.rec_gn + " : " + d.sig_gn;
        });
      // discard removed nodes
      info_box.exit().remove();
      // print the data again in case it changed
      info_box.text(function (d) {
        return d.rec_gn + " : " + d.sig_gn;
      });
    },

    conSeqUpdate(row_height) {
      var cScale = signprotmat.d3.cScale();
      var seqsigTip = d3
        .tip()
        .attr("class", "d3-tip")
        .html(function (d) {
          return (
            "Generic Residue No.: " +
            d.gn +
            "<br>" +
            "Feature: " +
            d.feature +
            "<br>" +
            "Length: " +
            d.length +
            "<br>" +
            // "Score: " +
            // d.expl +
            // "<br>" +
            "Frequency: " +
            d.freq +
            "<br>"
          );
        });
      // create selection and bind data
      var con_seq_mat = d3
        .select("g#con_seq_mat")
        .selectAll("g")
        .data(Object.values(con_seq))
        .join("g");
      var svg = d3.select("svg.svg-content.seqsig");
      var t = svg.transition().duration(750);

      // adding the labels and collapse/hide buttons
      var label_area = d3
        .select("g#con_seq_mat")
        .append("g")
        .attr("id", "labels");

      var gap = 40;

      // Labels
      label_area
        .append("text")
        .attr("class", "con_seq_label")
        .attr("x", -2)
        .attr("y", 75 + row_height * 1.5 + gap)
        .attr("dy", row_height / 2)
        .text("Property Consensus");

      label_area
        .append("text")
        .attr("class", "con_seq_label")
        .attr("x", -2)
        .attr("y", 75 + row_height * 2.5 + gap)
        .attr("dy", row_height / 4)
        .text("Length");

      label_area
        .append("text")
        .attr("class", "con_seq_label")
        .attr("x", -2)
        .attr("y", 75 + row_height * 3 + gap)
        .attr("dy", row_height / 4)
        .text("Property Conservation");

      label_area
        .append("text")
        .attr("class", "con_seq_label")
        .attr("x", -2)
        .attr("y", 75 + row_height * 3.5 + gap)
        .attr("dy", row_height / 4)
        .text("Interaction Conservation");

      label_area
        .append("text")
        .attr("class", "con_seq_label")
        .attr("x", -2)
        .attr("y", 75)
        .attr("dy", row_height / 2)
        .text("Sequence Consensus");

      label_area
        .append("text")
        .attr("class", "con_seq_label")
        .attr("x", -2)
        .attr("y", 75 + row_height)
        .attr("dy", row_height / 4)
        .text("Sequence Conservation");

      label_area
        .append("text")
        .attr("class", "con_seq_label")
        .attr("x", -2)
        .attr("y", 75 + row_height * 1.5)
        .attr("dy", row_height / 4)
        .text("Interaction Conservation");

      // Expand Buttons
      var tmp = label_area
        .append("g")
        .attr("id", "prop_con_expand")
        .on("click", function () {
          signprotmat.d3.toggle_prop_count();
        });

      tmp
        .append("rect")
        .attr("class", "expand")
        .attr("x", -180)
        .attr("y", 75 + row_height * 1.5 + gap)
        .attr("rx", "4")
        .attr("width", "50px")
        .attr("height", row_height);

      tmp
        .append("text")
        .attr("class", "expand")
        .attr("x", -180)
        .attr("y", 75 + row_height * 1.5)
        .attr("dx", "25px")
        .attr("dy", row_height / 2 + 3 + gap)
        .text("+");

      // tmp = label_area
      //     .append('g')
      //     .attr('id', 'seq_con_expand')
      //     .on('click', function(){signprotmat.d3.toggle_acid_count()})
      //
      // tmp.append('rect')
      //     .attr('class', 'expand')
      //     .attr('x', -180)
      //     .attr('y', 75 + row_height * 2)
      //     .attr('rx', '4')
      //     .attr('width', '50px')
      //     .attr('height', row_height)
      //
      // tmp.append('text')
      //     .attr('class', 'expand')
      //     .attr('x', -180)
      //     .attr('y', 75 + row_height * 2)
      //     .attr('dx', '25px')
      //     .attr("dy", row_height / 2 + 3)
      //     .text('+')

      // PROPERTY CODES
      con_seq_mat
        .selectAll("rect.res_rect")
        .data(function (d) {
          return d;
        })
        .join(
          function (enter) {
            return (
              enter
                .append("rect")
                .attr("class", "res_rect")
                .call(seqsigTip)
                .on("mouseover", function (d) {
                  seqsigTip.show(d);
                })
                .on("mouseout", function (d) {
                  seqsigTip.hide();
                })
                .style("fill", function (d) {
                  var gcol = signprotmat.d3.fScaleColor(d.feature_code);
                  if (typeof gcol !== "undefined") {
                    return gcol.bg_color;
                  } else {
                    return null;
                  }
                })
                //.style("fill", '#ffffff')
                //.call(enter => enter.transition(t)
                //.style("fill", function(d: any) {
                //const gcol = signprotmat.d3.fScaleColor(d.feature_code);
                //if (typeof gcol !== "undefined") {
                //return gcol.bg_color;
                //} else {
                //return null;
                //}
                //}))
                .style("stroke", "black")
                .attr("x", function (d) {
                  return xScale(d.gn) - xScale.step() / 2;
                })
                .attr("y", 75 + row_height * 1.5 + gap)
                .attr("width", xScale.step())
                .attr("height", row_height)
            );
          },
          function (update) {
            return (
              update
                .attr("x", function (d) {
                  return xScale(d.gn) - xScale.step() / 2;
                })
                .attr("y", 75 + row_height * 1.5 + gap)
                //.style("fill", '#ffffff')
                .call(function (update) {
                  return update.transition(t).style("fill", function (d) {
                    var gcol = signprotmat.d3.fScaleColor(d.feature_code);
                    if (typeof gcol !== "undefined") {
                      return gcol.bg_color;
                    } else {
                      return null;
                    }
                  });
                })
            );
          },
          function (exit) {
            return exit.remove();
          }
        );

      con_seq_mat
        .selectAll("text.res_label")
        .data(function (d) {
          return d;
        })
        .join(
          function (enter) {
            return enter
              .append("text")
              .attr("class", "res_label")
              .attr("text-anchor", "middle")
              .attr("x", function (d) {
                return xScale(d.gn);
              })
              .attr("y", 75 + row_height * 1.5 + gap)
              .attr("dy", row_height / 2)
              .style("fill", function (d) {
                var gcol = signprotmat.d3.fScaleColor(d.feature_code);
                if (typeof gcol !== "undefined") {
                  if (typeof gcol.font_color !== "undefined") {
                    return gcol.font_color;
                  } else {
                    return "#000000";
                  }
                } else {
                  return "#000000";
                }
              })
              .text(function (d) {
                return d.feature_code;
              });
          },
          function (update) {
            return update
              .attr("x", function (d) {
                return xScale(d.gn);
              })
              .attr("y", 75 + row_height * 1.5 + gap)
              .style("fill", function (d) {
                var gcol = signprotmat.d3.fScaleColor(d.feature_code);
                if (typeof gcol !== "undefined") {
                  if (typeof gcol.font_color !== "undefined") {
                    return gcol.font_color;
                  } else {
                    return "#000000";
                  }
                } else {
                  return "#000000";
                }
              })
              .text(function (d) {
                return d.feature_code;
              });
          },
          function (exit) {
            return exit.remove();
          }
        );

      // PROPERTY LENGTH
      con_seq_mat
        .selectAll("rect.len_rect")
        .data(function (d) {
          return d;
        })
        .join(
          function (enter) {
            return enter
              .append("rect")
              .attr("class", "len_rect")
              .style("fill", "#ffffff")
              .style("stroke", "black")
              .attr("x", function (d) {
                return xScale(d.gn) - xScale.step() / 2;
              })
              .attr("y", 75 + row_height * 2.5 + gap)
              .attr("width", xScale.step())
              .attr("height", row_height / 2);
          },
          function (update) {
            return update
              .attr("x", function (d) {
                return xScale(d.gn) - xScale.step() / 2;
              })
              .attr("y", 75 + row_height * 2.5 + gap);
          },
          function (exit) {
            return exit.remove();
          }
        );

      con_seq_mat
        .selectAll("text.len_label")
        .data(function (d) {
          return d;
        })
        .join(
          function (enter) {
            return enter
              .append("text")
              .attr("class", "len_label")
              .attr("text-anchor", "middle")
              .attr("x", function (d) {
                return xScale(d.gn);
              })
              .attr("y", 75 + row_height * 2.5 + gap)
              .attr("dy", row_height / 4)
              .style("fill", "#000000")
              .text(function (d) {
                return d.length;
              });
          },
          function (update) {
            return update
              .attr("x", function (d) {
                return xScale(d.gn);
              })
              .attr("y", 75 + row_height * 2.5 + gap)
              .style("fill", "#000000")
              .text(function (d) {
                return d.length;
              });
          },
          function (exit) {
            return exit.remove();
          }
        );

      // Property Sequence CONSERVATION
      con_seq_mat
        .selectAll("rect.cons_rect")
        .data(function (d) {
          return d;
        })
        .join(
          function (enter) {
            return enter
              .append("rect")
              .attr("class", "cons_rect")
              .style("fill", function (d) {
                if (d.cons === -1) {
                  return "#ffffff";
                } else {
                  return cScale(d.freq);
                }
              })
              .style("stroke", "black")
              .attr("x", function (d) {
                return xScale(d.gn) - xScale.step() / 2;
              })
              .attr("y", 75 + row_height * 3 + gap)
              .attr("width", xScale.step())
              .attr("height", row_height / 2);
          },
          function (update) {
            return update
              .style("fill", function (d) {
                if (d.cons === -1) {
                  return "#ffffff";
                } else {
                  return cScale(d.freq);
                }
              })
              .attr("x", function (d) {
                return xScale(d.gn) - xScale.step() / 2;
              })
              .attr("y", 75 + row_height * 3 + gap);
          },
          function (exit) {
            return exit.remove();
          }
        );

      con_seq_mat
        .selectAll("text.cons_label")
        .data(function (d) {
          return d;
        })
        .join(
          function (enter) {
            return enter
              .append("text")
              .attr("class", "cons_label")
              .attr("text-anchor", "middle")
              .attr("x", function (d) {
                return xScale(d.gn);
              })
              .attr("y", 75 + row_height * 3 + gap)
              .attr("dy", row_height / 4)
              .style("fill", function (d) {
                if (Math.abs(d.freq) >= 50) {
                  return "#eaeaea";
                } else if (Math.abs(d.freq) < 50) {
                  return "#000000";
                }
              })
              .text(function (d) {
                return d.freq;
              });
          },
          function (update) {
            return update
              .attr("x", function (d) {
                return xScale(d.gn);
              })
              .attr("y", 75 + row_height * 3 + gap)
              .style("fill", function (d) {
                if (Math.abs(d.freq) >= 50) {
                  return "#eaeaea";
                } else if (Math.abs(d.freq) < 50) {
                  return "#000000";
                }
              })
              .text(function (d) {
                return d.freq;
              });
          },
          function (exit) {
            return exit.remove();
          }
        );

      // Interaction conservation
      con_seq_mat
        .selectAll("rect.prop_seq.cons_rect")
        .data(function (d) {
          return d;
        })
        .join(
          function (enter) {
            return enter
              .append("rect")
              .attr("class", "prop_seq cons_rect")
              .style("fill", function (d) {
                return cScale(d.int_freq);
              })
              .style("stroke", "black")
              .attr("x", function (d) {
                return xScale(d.gn) - xScale.step() / 2;
              })
              .attr("y", 75 + row_height * 3.5 + gap)
              .attr("width", xScale.step())
              .attr("height", row_height / 2);
          },
          function (update) {
            return update
              .style("fill", function (d) {
                return cScale(d.int_freq);
              })
              .attr("x", function (d) {
                return xScale(d.gn) - xScale.step() / 2;
              })
              .attr("y", 75 + row_height * 3.5 + gap);
          },
          function (exit) {
            return exit.remove();
          }
        );

      con_seq_mat
        .selectAll("text.prop_seq.cons_label")
        .data(function (d) {
          return d;
        })
        .join(
          function (enter) {
            return enter
              .append("text")
              .attr("class", "prop_seq cons_label")
              .attr("text-anchor", "middle")
              .attr("x", function (d) {
                return xScale(d.gn);
              })
              .attr("y", 75 + row_height * 3.5 + gap)
              .attr("dy", row_height / 4)
              .style("fill", function (d) {
                if (d.int_freq >= 50) {
                  return "#eaeaea";
                } else {
                  return "#000000";
                }
              })
              .text(function (d) {
                return d.int_freq;
              });
          },
          function (update) {
            return update
              .attr("x", function (d) {
                return xScale(d.gn);
              })
              .attr("y", 75 + row_height * 3.5 + gap)
              .style("fill", function (d) {
                if (d.int_freq >= 50) {
                  return "#eaeaea";
                } else {
                  return "#000000";
                }
              })
              .text(function (d) {
                return d.int_freq;
              });
          },
          function (exit) {
            return exit.remove();
          }
        );

      // AMINO ACIDS
      con_seq_mat
        .selectAll("rect.amino_rect")
        .data(function (d) {
          return d;
        })
        .join(
          function (enter) {
            return enter
              .append("rect")
              .attr("class", "amino_rect")
              .style("fill", function (d) {
                return signprotmat.d3.resScaleColor(d.aa)["bg_color"];
              })
              .style("stroke", "black")
              .attr("x", function (d) {
                return xScale(d.gn) - xScale.step() / 2;
              })
              .attr("y", 75)
              .attr("width", xScale.step())
              .attr("height", row_height);
          },
          function (update) {
            return update
              .attr("x", function (d) {
                return xScale(d.gn) - xScale.step() / 2;
              })
              .attr("y", 75)
              .call(function (update) {
                return update.transition(t).style("fill", function (d) {
                  return signprotmat.d3.resScaleColor(d.aa)["bg_color"];
                });
              });
          },
          function (exit) {
            return exit.remove();
          }
        );

      con_seq_mat
        .selectAll("text.amino_label")
        .data(function (d) {
          return d;
        })
        .join(
          function (enter) {
            return enter
              .append("text")
              .attr("class", "amino_label")
              .attr("text-anchor", "middle")
              .attr("x", function (d) {
                return xScale(d.gn);
              })
              .attr("y", 75)
              .attr("dy", row_height / 2)
              .style("fill", function (d) {
                return signprotmat.d3.resScaleColor(d.aa)["fg_color"];
              })
              .text(function (d) {
                return d.aa;
              });
          },
          function (update) {
            return update
              .attr("x", function (d) {
                return xScale(d.gn);
              })
              .attr("y", 75)
              .style("fill", function (d) {
                return signprotmat.d3.resScaleColor(d.aa)["fg_color"];
              })
              .text(function (d) {
                return d.aa;
              });
          },
          function (exit) {
            return exit.remove();
          }
        );

      // AMINO CONSERVATION
      con_seq_mat
        .selectAll("rect.amino.cons_rect")
        .data(function (d) {
          return d;
        })
        .join(
          function (enter) {
            return enter
              .append("rect")
              .attr("class", "amino cons_rect")
              .style("fill", function (d) {
                return cScale(d.aa_cons);
              })
              .style("stroke", "black")
              .attr("x", function (d) {
                return xScale(d.gn) - xScale.step() / 2;
              })
              .attr("y", 75 + row_height)
              .attr("width", xScale.step())
              .attr("height", row_height / 2);
          },
          function (update) {
            return update
              .attr("x", function (d) {
                return xScale(d.gn) - xScale.step() / 2;
              })
              .attr("y", 75 + row_height)
              .call(function (update) {
                return update.transition(t).style("fill", function (d) {
                  return cScale(d.aa_cons);
                });
              });
          },
          function (exit) {
            return exit.remove();
          }
        );

      con_seq_mat
        .selectAll("text.amino.cons_label")
        .data(function (d) {
          return d;
        })
        .join(
          function (enter) {
            return enter
              .append("text")
              .attr("class", "amino cons_label")
              .attr("text-anchor", "middle")
              .attr("x", function (d) {
                return xScale(d.gn);
              })
              .attr("y", 75 + row_height)
              .attr("dy", row_height / 4)
              .style("fill", function (d) {
                if (Math.abs(d.aa_cons) >= 50) {
                  return "#eaeaea";
                } else if (Math.abs(d.aa_cons) < 50) {
                  return "#000000";
                }
              })
              .text(function (d) {
                return d.aa_cons;
              });
          },
          function (update) {
            return update
              .attr("x", function (d) {
                return xScale(d.gn);
              })
              .attr("y", 75 + row_height)
              .style("fill", function (d) {
                if (Math.abs(d.aa_cons) >= 50) {
                  return "#eaeaea";
                } else if (Math.abs(d.aa_cons) < 50) {
                  return "#000000";
                }
              })
              .text(function (d) {
                return d.aa_cons;
              });
          },
          function (exit) {
            return exit.remove();
          }
        );

      // Interaction conservation
      con_seq_mat
        .selectAll("rect.seq.cons_rect")
        .data(function (d) {
          return d;
        })
        .join(
          function (enter) {
            return enter
              .append("rect")
              .attr("class", "seq cons_rect")
              .style("fill", function (d) {
                return cScale(d.int_freq);
              })
              .style("stroke", "black")
              .attr("x", function (d) {
                return xScale(d.gn) - xScale.step() / 2;
              })
              .attr("y", 75 + row_height * 1.5)
              .attr("width", xScale.step())
              .attr("height", row_height / 2);
          },
          function (update) {
            return update
              .style("fill", function (d) {
                return cScale(d.int_freq);
              })
              .attr("x", function (d) {
                return xScale(d.gn) - xScale.step() / 2;
              })
              .attr("y", 75 + row_height * 1.5);
          },
          function (exit) {
            return exit.remove();
          }
        );

      con_seq_mat
        .selectAll("text.seq.cons_rect")
        .data(function (d) {
          return d;
        })
        .join(
          function (enter) {
            return enter
              .append("text")
              .attr("class", "seq cons_label")
              .attr("text-anchor", "middle")
              .attr("x", function (d) {
                return xScale(d.gn);
              })
              .attr("y", 75 + row_height * 1.5)
              .attr("dy", row_height / 4)
              .style("fill", function (d) {
                if (d.int_freq >= 50) {
                  return "#eaeaea";
                } else {
                  return "#000000";
                }
              })
              .text(function (d) {
                return d.int_freq;
              });
          },
          function (update) {
            return update
              .attr("x", function (d) {
                return xScale(d.gn);
              })
              .attr("y", 75 + row_height * 1.5)
              .style("fill", function (d) {
                if (d.int_freq >= 50) {
                  return "#eaeaea";
                } else {
                  return "#000000";
                }
              })
              .text(function (d) {
                return d.int_freq;
              });
          },
          function (exit) {
            return exit.remove();
          }
        );
    },

    draw_seq_sig(data_in, svg, xScale) {
      var data = data_in.feat;

      var cScale = signprotmat.d3.cScale();
      var feats = [];

      for (var _i = 0, _a = Object.keys(data); _i < _a.length; _i++) {
        var key = _a[_i];
        if (xScale(key) == null) {
          data = _.omit(data, key);
        }
      }

      var col_lengths = [];
      for (var _b = 0, _c = Object.keys(data); _b < _c.length; _b++) {
        var elem = _c[_b];
        col_lengths.push(data[elem].length);
      }

      var row_height = 30;
      var area_height = _.max(col_lengths) * row_height;
      var seqsigTip = d3
        .tip()
        .attr("class", "d3-tip")
        .html(function (d) {
          return (
            "Generic Residue No.: " +
            d.gn +
            "<br>" +
            "Feature: " +
            d.feature +
            "<br>" +
            "Length: " +
            d.length +
            "<br>" +
            // "Score: " +
            // d.expl +
            // "<br>" +
            "Frequency: " +
            d.freq +
            "<br>"
          );
        });

      var data = data_in.feat_ungrouped;
      var fScale = signprotmat.d3.fScale(data);

      data.forEach((d) => {
        const length_text = d.length != "" ? " (" + d.length + ")" : "";
        feats.push({
          code: d.feature_code,
          feature: d.feature,
          length: d.length,
          comb: d.feature + length_text,
        });
      });
      let uniq_feats = _.uniqBy(feats, "comb");

      // filter out NA generic numbers based on xScale
      data = _.filter(data, function (d) {
        return xScale(d.gn);
      });

      let row = svg
        .append("g")
        .attr("id", "seqsig_group")
        .attr("class", "collapse")
        .append("g")
        .attr("id", "seqsig_feature")
        .attr("transform", "translate(" + 0 + "," + 140 + ")")
        .selectAll("text")
        .data(uniq_feats)
        .enter();

      //feature row labels
      row
        .append("text")
        .attr("class", "y seq_label")
        .attr("x", -30 - xScale.step())
        .attr("y", function (d) {
          return fScale(d.comb) - fScale.step() / 2;
        })
        .attr("text-anchor", "end")
        .attr("dy", 75)
        .text(function (d) {
          return d.feature;
        });

      //feature length row labels
      row
        .append("text")
        .attr("class", "y seq_label")
        .attr("x", -10 - xScale.step())
        .attr("y", function (d) {
          return fScale(d.comb) - fScale.step() / 2;
        })
        .attr("text-anchor", "end")
        .attr("dy", 75)
        .text(function (d) {
          return d.length;
        });

      // feature color code rectangles
      row
        .append("rect")
        .style("fill", function (d) {
          const gcol = signprotmat.d3.fScaleColor(d.code);
          if (typeof gcol !== "undefined") {
            return gcol.bg_color;
          } else {
            return null;
          }
        })
        .style("stroke", "black")
        .attr("x", -xScale.step())
        .attr("y", function (d) {
          return 75 + fScale(d.comb) - fScale.step();
        })
        .attr("width", xScale.step())
        .attr("height", fScale.step());

      // feature code label text
      row
        .append("text")
        .attr("class", "y seq_label")
        .attr("text-anchor", "middle")
        .attr("x", -xScale.step() / 2)
        .attr("y", function (d) {
          return 75 + fScale(d.comb) - fScale.step() / 2;
        })
        .style("fill", (d) => {
          const gcol = signprotmat.d3.fScaleColor(d.code);
          if (typeof gcol !== "undefined") {
            if (typeof gcol.font_color !== "undefined") {
              return gcol.font_color;
            } else {
              return "#000000";
            }
          } else {
            return "#000000";
          }
        })
        .text(function (d) {
          return d.code;
        });

      // generating the white backdrop for all the properties
      svg
        .select("g#seqsig_group")
        .append("g")
        .attr("id", "seqsig_mat")
        .attr("transform", "translate(" + -xScale.step() / 2 + "," + 140 + ")")
        .append("rect")
        .attr("class", "border-bg")
        .style("fill", "#ffffff")
        .attr("x", xScale.step() / 2)
        .attr("y", 75)
        .attr("width", xScale.range()[1] - xScale.step())
        .attr("height", fScale.range()[1] - fScale.step());

      let each_res = svg
        .select("g#seqsig_mat")
        .selectAll("text")
        .data(data)
        .enter()
        .append("g")
        .call(seqsigTip)
        .on("mouseover", function (d) {
          if (d.freq !== 0) {
            seqsigTip.show(d);
          }
        })
        .on("mouseout", function (d) {
          seqsigTip.hide();
        });

      // the rectangles, colored by conservation
      each_res
        .append("rect")
        .attr("class", "res_rect")
        .style("fill", function (d) {
          if (d.cons <= 0) {
            return "none";
          } else {
            return cScale(d.freq);
          }
        })
        .attr("x", (d) => xScale(d.gn) - xScale.step() / 2)
        .attr("y", function (d) {
          const length_text = d.length != "" ? " (" + d.length + ")" : "";
          const comb = d.feature + length_text;
          return 75 + fScale(comb) - fScale.step();
        })
        .attr("width", xScale.step())
        .attr("height", fScale.step());

      // adding the frequency text to each rectangle
      each_res
        .append("text")
        .attr("class", "res_label")
        .attr("x", function (d) {
          return xScale(d.gn);
        })
        .attr("y", function (d) {
          const length_text = d.length != "" ? " (" + d.length + ")" : "";
          const comb = d.feature + length_text;
          return fScale(comb) - fScale.step() / 2;
        })
        .style("fill", function (d) {
          if (Math.abs(d.freq) >= 50) {
            return "#eaeaea";
          } else if (Math.abs(d.freq) < 50) {
            return "#000000";
          }
        })
        .attr("text-anchor", "middle")
        .attr("dy", 75)
        .text(function (d) {
          return d.freq;
        });
      // .text(function(d) { return _.round(d.freq/100, 1) });

      // Setting up the consensus area
      d3.select("svg.svg-content.seqsig")
        .select("g")
        .append("g")
        .attr("id", "con_seq_mat")
        .attr("transform", "translate(" + -xScale.step() / 2 + "," + -70 + ")");

      // adding the gn labels
      d3.select("svg.svg-content.seqsig")
        .select("g")
        .append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(" + -xScale.step() / 2 + "," + 170 + ")")
        .call(xAxis)
        .selectAll("text")
        .attr("text-anchor", "end")
        .attr("font-size", "12px")
        .attr("dx", "-5px")
        .attr("dy", "-5px")
        .attr("transform", "rotate(-90)");

      // removing the axis line
      svg.select(".x.axis").selectAll("path").remove();

      // putting a black border around the signature
      d3.select("g#seqsig_mat")
        .append("rect")
        .attr("class", "border")
        .style("stroke", "black")
        .style("fill", "none")
        .attr("x", xScale.step() / 2)
        .attr("y", 75)
        .attr("width", xScale.range()[1] - xScale.step())
        .attr("height", fScale.range()[1] - fScale.step());

      // reducing the viewbox height to the required height
      var viewbox_svg = d3.select(".svg-content.seqsig");
      var viewbox = viewbox_svg.attr("viewBox");
      var viewbox_1 = viewbox.slice(0, 9);
      // var new_vb = viewbox_1 + (area_height + Math.round(area_height / 4) + 150);
      var new_vb = viewbox_1 + 210;
      viewbox_svg.attr("viewBox", new_vb);
    },

    toggle_prop_count() {
      // check to see if the table is visible
      bool_visible = d3.select("g#seqsig_group").classed("in");

      // toggle the table visiblity state
      d3.select("g#seqsig_group").classed("in", !bool_visible);

      //set viewbox to correct height
      var viewbox_svg = d3.select(".svg-content.seqsig");
      var viewbox = viewbox_svg.attr("viewBox");
      var viewbox_1 = viewbox.slice(0, 9);
      var new_vb = viewbox_1 + (!bool_visible ? 1510 : 210);
      viewbox_svg.attr("viewBox", new_vb);

      // redraw the container to fix viewbox display height
      d3.select("div#seqsig-svg").classed("in", false);
      setTimeout(() => d3.select("div#seqsig-svg").classed("in", true), 50);
    },

    toggle_acid_count() {
      // check to see if the table is visible
      bool_visible = d3.select("g#acid_group").classed("in");

      // toggle the table visiblity state
      d3.select("g#acid_group").classed("in", !bool_visible);

      //set viewbox to correct height
      var viewbox_svg = d3.select(".svg-content.seqsig");
      var viewbox = viewbox_svg.attr("viewBox");
      var viewbox_1 = viewbox.slice(0, 9);
      var new_vb = viewbox_1 + (!bool_visible ? 1510 : 210);
      viewbox_svg.attr("viewBox", new_vb);

      // redraw the container to fix viewbox display height
      d3.select("div#seqsig-svg").classed("in", false);
      setTimeout(() => d3.select("div#seqsig-svg").classed("in", true), 10);
    },

    draw_seq_cons(data_in, svg, xScale, xAxis, sigmatch) {
      var data = data_in.cons;
      var fScale = signprotmat.d3.fScale(data);
      var cScale = signprotmat.d3.cScale(data);
      var uniq_feats = _.uniq(_.map(data, "feature"));
      // filter out NA generic numbers based on xScale
      data = _.filter(data, function (d) {
        return xScale(d.gn);
      });
      var conseqTip = d3
        .tip()
        .attr("class", "d3-tip")
        .html(function (d) {
          return (
            "Generic Residue No.: " +
            d.gn +
            "<br>" +
            "Feature: " +
            d.feature +
            "<br>" +
            "Length: " +
            d.length +
            "<br>" +
            "Score: " +
            d.score +
            "<br>"
          );
        });
      if (sigmatch) {
        svg = svg
          .append("g")
          .attr("id", "sigmatch_frame")
          .attr(
            "transform",
            "translate(" + margin.left + "," + (margin.top + 140) + ")"
          );
        svg
          .append("g")
          .attr("id", "sigmatch_mat")
          .attr(
            "transform",
            "translate(" + -xScale.step() / 2 + "," + -30 + ")"
          )
          .append("rect")
          .attr("class", "border-bg")
          .style("fill", "#ffffff")
          .attr("x", xScale.step() / 2)
          .attr("y", 75)
          .attr("width", xScale.range()[1] - xScale.step())
          .attr("height", 75);
        svg
          .append("text")
          .attr("class", "y seq_label")
          .attr("text-anchor", "end")
          .attr("x", -10)
          .attr("y", 65)
          .text("Property");
        svg
          .append("text")
          .attr("class", "y seq_label")
          .attr("text-anchor", "end")
          .attr("x", -10)
          .attr("y", 102)
          .text("Conservation");
        var group = "g#sigmatch_mat";
      } else {
        svg
          .append("g")
          .attr("id", "conseq_mat")
          .attr(
            "transform",
            "translate(" + -xScale.step() / 2 + "," + -30 + ")"
          )
          .append("rect")
          .attr("class", "border-bg")
          .style("fill", "#ffffff")
          .attr("x", xScale.step() / 2)
          .attr("y", 75)
          .attr("width", xScale.range()[1] - xScale.step())
          .attr("height", 75);
        svg
          .append("g")
          .attr("class", "x axis")
          .attr("transform", "translate(" + -xScale.step() / 2 + "," + 35 + ")")
          .call(xAxis)
          .selectAll("text")
          .attr("text-anchor", "start")
          .attr("font-size", "12px")
          .attr("dx", "-5px")
          .attr("dy", "-5px")
          .attr("transform", "rotate(-90)");
        svg.select(".x.axis").selectAll("path").remove();
        svg
          .append("text")
          .attr("class", "y seq_label")
          .attr("text-anchor", "end")
          .attr("x", -10)
          .attr("y", 65)
          .text("Property");
        svg
          .append("text")
          .attr("class", "y seq_label")
          .attr("text-anchor", "end")
          .attr("x", -10)
          .attr("y", 102)
          .text("Conservation");
        var group = "g#conseq_mat";
      }

      var each_res = svg
        .select(group)
        .selectAll("text")
        .data(data)
        .enter()
        .append("g")
        .call(conseqTip)
        .on("mouseover", function (d) {
          conseqTip.show(d);
        })
        .on("mouseout", function (d) {
          conseqTip.hide();
        });

      // the rectangles, colored by feature
      each_res
        .append("rect")
        .attr("class", "res_rect")
        .style("fill", function (d) {
          var gcol = signprotmat.d3.fScaleColor(d.code);
          if (typeof gcol !== "undefined") {
            return gcol.bg_color;
          } else {
            return null;
          }
        })
        .attr("x", function (d) {
          return xScale(d.gn) - xScale.step() / 2;
        })
        .attr("y", function (d) {
          return 75;
        })
        .attr("width", xScale.step())
        .attr("height", 37.5);

      // the rectangles, colored by conservation
      each_res
        .append("rect")
        .attr("class", "res_rect")
        .style("fill", function (d) {
          if (d.cons === -1) {
            return "#ffffff";
          } else {
            return cScale(d.score);
          }
        })
        .attr("x", function (d) {
          return xScale(d.gn) - xScale.step() / 2;
        })
        .attr("y", function (d) {
          return 75 + 37.5;
        })
        .attr("width", xScale.step())
        .attr("height", 37.5);

      // adding the feature text to each rectangle
      each_res
        .append("text")
        .attr("class", "res_label")
        // .attr("x", (d: any) => xScale(d.gn))
        // .attr("y", (d: any) => 50)
        .attr(
          "transform",
          function (d) {
            return "translate(" + xScale(d.gn) + ",93.75)";
          } // + "rotate(270)"
        )
        .style("fill", function (d) {
          var gcol = signprotmat.d3.fScaleColor(d.code);
          if (typeof gcol !== "undefined") {
            if (typeof gcol.font_color !== "undefined") {
              return gcol.font_color;
            } else {
              return "#000000";
            }
          } else {
            return "#000000";
          }
        })
        .attr("text-anchor", "middle")
        .text(function (d) {
          return d.code;
        });

      // adding the conservation value to each rectangle
      each_res
        .append("text")
        .attr("class", "res_label")
        // .attr("x", (d: any) => xScale(d.gn))
        // .attr("y", (d: any) => 50)
        .attr(
          "transform",
          function (d) {
            return (
              "translate(" + xScale(d.gn) + "," + (75 + 37.5 + 37.5 / 2) + ")"
            );
          } // + "rotate(270)"
        )
        .style("fill", function (d) {
          if (Math.abs(d.score) >= 50) {
            return "#eaeaea";
          } else if (Math.abs(d.score) < 50) {
            return "#000000";
          }
        })
        .attr("text-anchor", "middle")
        .text(function (d) {
          return d.score;
        });

      // putting a black border around the signature
      d3.select(group)
        .append("rect")
        .attr("class", "border")
        .style("stroke", "black")
        .style("fill", "none")
        .attr("x", xScale.step() / 2)
        .attr("y", 75)
        .attr("width", xScale.range()[1] - xScale.step())
        .attr("height", 75);
    },
  },
};
