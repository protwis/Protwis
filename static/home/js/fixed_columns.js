/*eslint complexity: ["error", 20]*/
var oTable = [];
function update_text_in_modal() {

    var pdbs = [];
    var mode = $("ul#mode_nav").find("li.active").find("a").text().trim();
    var ModalpdbsCountSelector;
    var ExternalModalpdbCounter;
    var total_selected;
    var pdbsInputSelector;
    var pdbsCountSelector;
    var selected_visible;
    var selector;
    var group;
    group = $(".tableview:visible").attr("group-number");
    if (group) mode = mode + group;

    $(".pdb_selected", oTable[mode].cells().nodes()).each(function() {
        if ($(this).prop("checked")) {
            $(this).closest("tr").addClass("selected");
            $(this).closest(".dataTables_scroll").find("#overlay_" + $(this).attr("id")).addClass("selected");
            selector = $(this).closest(".dataTables_scrollBody").find("#overlay_" + $(this).attr("id")).find("checkbox");
            selector.prop("checked", !selector.prop("checked"));
            pdbs.push($(this).attr("id"));
        } else {
            $(this).closest("tr").removeClass("selected");
            $(this).closest(".dataTables_scroll").find("#overlay_" + $(this).attr("id")).removeClass("selected");
        }
    });

    mode = $("ul#mode_nav").find("li.active").find("a").text().trim();
    if (mode === "Single set of structures" || $("#single-group-tree-tab").length) {
        total_selected = pdbs.length;
        selected_visible = $(".dataTables_scrollBody:visible .pdb_selected:checked").length;
        ModalpdbsCountSelector = "#single-crystal-group-pdbs-modal-text";
        ExternalModalpdbCounter = "#single-crystal-group-pdbs-modal-external-text";

        if (total_selected == selected_visible) {
            $(ModalpdbsCountSelector).html(total_selected + " structure(s) selected (ticked leftmost)");
            $(ExternalModalpdbCounter).html(total_selected + " structure(s) selected");
        } else {
            $(ModalpdbsCountSelector).html(total_selected + " structure(s) selected (ticket leftmost — " + (total_selected - selected_visible) + " currently filtered)");
            $(ExternalModalpdbCounter).html(total_selected + " structure(s) selected (" + (total_selected - selected_visible) + " currently filtered)");
        }

        pdbsInputSelector = "#single-crystal-group-tab .crystal-pdb";
        pdbsCountSelector = "#single-crystal-group-tab .crystal-count";
        $(pdbsInputSelector).val(JSON.stringify(pdbs));
        $(pdbsCountSelector).html(pdbs.length);

    } else if (mode === "Two sets of structures") {
        total_selected = pdbs.length;
        selected_visible = $(".pdb_selected:checked:visible").length;
        ModalpdbsCountSelector = "#two-crystal-group-pdbs-modal-" + group + "-text";

        if (total_selected == selected_visible) {
            $(ModalpdbsCountSelector).html(total_selected + " structure(s) selected");
        } else {
            $(ModalpdbsCountSelector).html(total_selected + " structure(s) selected (" + (total_selected - selected_visible) + " currently filtered)");
        }

        pdbsInputSelector = "#two-crystal-groups-tab .crystal-group-" + group + "-pdbs";
        pdbsCountSelector = "#two-crystal-groups-tab .crystal-count-" + group;
        $(pdbsInputSelector).val(JSON.stringify(pdbs));
        $(pdbsCountSelector).html(pdbs.length);

    }
}

function thisPDB(elem) {

    loading_anim = $('#interface-modal-table' + " .loading_overlay");
    loading_anim.show();
    setTimeout(() => {

        var mode = $("ul#mode_nav").find("li.active").find("a").text().trim();
        var ReceptorName = $(elem).attr("long");

        var referenceObject = $(elem).closest(".dataTables_scroll")

        // Fixed/frozen columns check
        if (elem.id.startsWith("overlaycheck_")){
            elem = $(".dataTables_scrollBody:visible").find("#" + elem.id.replace("overlaycheck_", ""))[0]

        $(elem).prop("checked", !elem.checked)
        } else {
        other = referenceObject.find("#overlaycheck_" + elem.id)
        other.prop("checked", elem.checked)
        }

        var pdbName = $(elem).attr("id");
        if (mode === "Single structure") {
            if ($("input", oTable[mode].cells().nodes()).filter(":checkbox").not(elem).length > 0) {
            $("input", oTable[mode].cells().nodes()).filter(":checkbox").not(elem).each(function(i,e)
                { if (referenceObject.find("#overlaycheck_" + e.id).length > 0)
                { referenceObject.find("#overlaycheck_" + e.id)[ 0 ].checked = false;
                }
                }
            ); // deselect from overlay
            $("input", oTable[mode].cells().nodes()).filter(":checkbox").not(elem).prop("checked", false); // deselect from original table
            }
            var pdbs = [];
            if ($(elem).prop("checked")) {
                pdbs.push(pdbName);
                // Update view
                $(".crystal-count:visible").html(ReceptorName + " - " + pdbName + " selected.");
            } else {
                // Update view
                $(".crystal-count:visible").html("No structure selected.");
            }
            $(".crystal-count:visible").parent().parent().find(".crystal-pdb").val(JSON.stringify(pdbs));
        }

        update_text_in_modal();
        loading_anim.hide();
    }, 100);
}


function resetselection(not_update = false, reset_filters = false) {
    var mode = $("ul#mode_nav").find("li.active").find("a").text().trim();
    var group = $(".tableview:visible").attr("group-number");
    if (group) mode = mode + group;

    $(".check_all:visible").prop("checked", false);

    $("input", oTable[mode].cells().nodes()).prop("checked", false);

    if (!not_update) update_text_in_modal();

    if (reset_filters) yadcf.exResetAllFilters(oTable[mode]);
}

function check_all(elem, button) {
    var show_all;
    var group;
    var mode = $("ul#mode_nav").find("li.active").find("a").text().trim();
    show_all = $(".check_all:visible").prop("checked");
    if (button) {
        if (show_all) {
            $(".check_all:visible").prop("checked", false);
            show_all = false;
        } else {
            $(".check_all:visible").prop("checked", true);
            show_all = true;
        }
    }

    if (mode === "Single set of structures" || $("#single-group-tree-tab").length) {
        var pdbs = [];

        // REMOVE EXISITING? Probably not, more logical that filtering adds more
        // $("input", oTable.cells().nodes()).prop("checked",false);

        if (show_all) {
            $(".pdb_selected:visible").prop("checked", true);
        } else {
            $(".pdb_selected:visible").prop("checked", false);
        }
    } else if (mode === "Two sets of structures") {
        group = $(elem).closest(".tableview").attr("group-number");
        if (group) mode = mode + group;
        var pdbs = [];

        if (show_all) {
            $(".pdb_selected:visible").prop("checked", true);
        } else {
            $(".pdb_selected:visible").prop("checked", false);
        }
    }
    if (button) {
        update_text_in_modal();
    }
}

function check_all_representatives() {
    var mode = $("ul#mode_nav").find("li.active").find("a").text().trim();
    var group = $(".tableview:visible").attr("group-number");
    if (group) mode = mode + group;
    $("input", oTable[mode].cells().nodes()).prop("checked", false);
    $('input[representative="Yes"]:visible').each(function() {
        $(this).prop("checked", true);
    });
    update_text_in_modal();
}

function check_all_distance_representatives() {
    var mode = $("ul#mode_nav").find("li.active").find("a").text().trim();
    var group = $(".tableview:visible").attr("group-number");
    if (group) mode = mode + group;
    $("input", oTable[mode].cells().nodes()).prop("checked", false);
    $('input[distance_representative="Yes"]:visible').each(function() {
        // pdbs.push($(this).attr("id"));
        $(this).prop("checked", true);
    });
    update_text_in_modal();
}

function check_all_class_representatives() {
    var mode = $("ul#mode_nav").find("li.active").find("a").text().trim();
    let group = $(".tableview:visible").attr("group-number");
    if (group) mode = mode + group;
    $("input", oTable[mode].cells().nodes()).prop("checked", false);
    $('input[class_consensus_based_representative="Yes"]:visible').each(function() {
        // pdbs.push($(this).attr('id'));
        $(this).prop("checked", true);
    });
    update_text_in_modal();
}

$.fn.dataTable.ext.order["dom-checkbox"] = function(settings, col) {
    return this.api().column(col, {
        order: "index"
    }).nodes().map(function(td, i) {
        return $("input", td).prop("checked") ? "1" : "0";
    });
};

function pastePDBs() {
    var mode = $("ul#mode_nav").find("li.active").find("a").text().trim();
    var group = $(".tableview:visible").attr("group-number");
    if (group) mode = mode + group;
    var pdbs = $(".pastePDBs:visible").val().toUpperCase().trim();
    pdbs = pdbs.split(/[ ,]+/);
    resetselection(1);
    $(".pdb_selected", oTable[mode].cells().nodes()).each(function() {
        var pdb = $(this).attr("id");
        check = $.inArray(pdb, pdbs);
        if (check !== -1) {
            $(this).prop("checked", true);
            pdbs.splice(check, 1);
        }
    });
    if (pdbs.length) {
        var popOverSettings = {
            placement: "bottom",
            container: "body",
            title: "One or more PDB(s) not found:", //Sepcify the selector here
            content: pdbs
        }
        $(".pastePDBs").popover(popOverSettings);
        $(".pastePDBs").popover("show");
        $(".pastePDBs").on("shown.bs.popover",function() {
         setTimeout(function() {
          $(".pastePDBs").popover("destroy");
          }, 3000);
        });
    } /*else {
        $(".pastePDBs").popover("destroy");
    }*/
    oTable[mode].order([
        [20, "desc"]
    ]);
    oTable[mode].draw();
}

function prepopulatePDBs() {
    var pdbsInputSelector;
    var pdb;
    var mode = $("ul#mode_nav").find("li.active").find("a").text().trim();
    var mode2 = $("ul#mode_nav").find("li.active").find("a").text().trim();
    var group = $(".tableview:visible").attr("group-number");
    if (group) mode = mode + group;
    if (mode2 === "Two sets of structures") {
        pdbsInputSelector = "#two-crystal-groups-tab .crystal-group-" + group + "-pdbs";
        var pdbs = JSON.parse($(pdbsInputSelector).val());
        $(".pdb_selected", oTable[mode].cells().nodes()).each(function () {
            pdb = $(this).attr("id");
            check = $.inArray(pdb, pdbs);
            if (check !== -1) {
                $(this).prop("checked", true);
                pdbs.splice(check, 1);
            }
        });

        oTable[mode].order([
            [0, "desc"]
        ]);
        oTable[mode].draw();
    } else if (mode2 === "Single set of structures") {
        pdbsInputSelector = "#single-crystal-group-tab .crystal-pdb";
        if ( $( pdbsInputSelector ).length ) {
          var pdbs = JSON.parse($(pdbsInputSelector).val());
          $(".pdb_selected", oTable[mode].cells().nodes()).each(function () {
              pdb = $(this).attr("id");
              check = $.inArray(pdb, pdbs);
              if (check !== -1) {
                  $(this).prop("checked", true);
                  pdbs.splice(check, 1);
              }
          });

          oTable[mode].order([
              [0, "desc"]
          ]);
          oTable[mode].draw();
        }
    }
}

function exportPDBs() {
    var mode = $("ul#mode_nav").find("li.active").find("a").text().trim();
    var group = $(".tableview:visible").attr("group-number");
    if (group) mode = mode + group;
    var pdbs = [];
    $(".pdb_selected:checked", oTable[mode].cells().nodes()).each(function() {
        pdbs.push($(this).attr("id"));
    });


    var textArea = document.createElement("textarea");
    textArea.value = pdbs;
    document.body.appendChild(textArea);
    textArea.focus();
    textArea.select();
    try {
        var successful = document.execCommand("copy");
        var msg = successful ? "Successful" : "Unsuccessful";
        $(".export_pdbs").html(msg);
        setTimeout("$('.export_pdbs').html('Export selected PDB codes');", 4000);
    } catch (err) {
        $(".export_pdbs").html("Oops, unable to copy");
    }

    document.body.removeChild(textArea);
}

function toggle_best(mode, index, value) {
    var filter_value;
    filter_value = value === "On" ? "Best" : "";
    yadcf.exFilterColumn(oTable[mode], [[index, filter_value]]);
}

function create_overlay(table_id) {
    // This function fires upon filtering, to update what rows to show as an overlay
    $(".overlay_table tbody tr").remove();
    var $target = $(".overlay_table tbody");
    $(table_id + " tbody tr").each(function() {
        var $tds = $(this).children(),
            $row = $("<tr id='overlay_" + $tds.eq(7).html() + "''></tr>");
        // $row.append($tds.eq(0).clone()).append($tds.eq(1).clone()).appendTo($target);
        $row.append($tds.eq(0).clone()).append($tds.eq(1).clone()).append($tds.eq(2).clone()).append($tds.eq(3).clone()).appendTo($target);
    });
    $(".overlay_table .border-right").removeClass("border-right");

    // rename checkboxes for overlay to link to the original and remove class
    $(".overlay_table tbody tr :checkbox").each(function(i, e){
      e.id = "overlaycheck_" + e.id;
      e.classList.remove("pdb_selected");
    });
    // rebind click event
    /*$(".structure_overlay tr").click(function(event) {
        console.log("clicking overlay tr");
        if (event.target.type !== "checkbox") {
            $(":checkbox", this).trigger("click");
            pdb_id = $(this).attr("id").split("_")[1];
            console.log(pdb_id);
            if ($(":checkbox", this).length === 0) {
                pdb_id = $(this).attr("id").split("_")[1];
                console.log("2", pdb_id, "#" + pdb_id + ":checkbox:visible", $("#" + pdb_id + ":checkbox:visible"));
                // $("#"+pdb_id+":checkbox:visible").find("tr").trigger("click");
            }
        }
    });*/
    $(".structure_overlay tr").click(function(event) {
        if (event.target.type !== "checkbox") {
            $(":checkbox", this).trigger("click");
            if ($(":checkbox", this).length === 0) {
                var pdb_id = $(this).attr("id").split("_")[1];
                var checkbox = $(".dataTables_scrollBody").find("#" + pdb_id);
                checkbox.trigger("click");
            }
        }
    });
    $(".structure_overlay tr").css("cursor", "pointer");
}

function showPDBtable(element) {
    var mode = $("ul#mode_nav").find("li.active").find("a").text().trim();
    var group = $(element + " .tableview").attr("group-number");
    if (group) mode = mode + group;
    // console.log(element,mode,group);
    var mode_without_space;
    mode_without_space = mode.replace(/ /g, "_");

    if (!$.fn.DataTable.isDataTable(element + " .tableview table")) {
        $(element + " #best_species").html('<div class2="pull-right">\
                                            <div class="btn-group btn-toggle" column="7" mode="'+ mode +'"> \
                                            <button class="btn btn-xs btn-default" value="On">&nbsp;</button> \
                                            <button class="btn btn-xs btn-primary active" value="Off">Off</button> \
                                            </div> \
                                            <span class="pull-right glyphicon glyphicon-info-sign" title="Human or closest species (for each receptor and state)" style="font-size: 15px;cursor: pointer;"></span> \
                                            </div>');

        $(element + " #best_res").html('<div class2="pull-right"> \
        <div class="btn-group btn-toggle" column="12" mode="'+ mode +'"> \
                                            <button class="btn btn-xs btn-default" value="On">&nbsp;</button> \
                                            <button class="btn btn-xs btn-primary active" value="Off">Off</button> \
                                            </div> \
                                            <span class="pull-right glyphicon glyphicon-info-sign" title="Highest resolution (for each receptor, species and state)" style="font-size: 15px;cursor: pointer;"></span> \
                                            </div>');


        $(element + " .btn-toggle").click(function() {
            $(this).find(".btn").toggleClass("active");
            $(this).find(".btn").toggleClass("btn-primary");
            $(this).find(".btn").toggleClass("btn-default");
            $(this).find(".active").html($(this).find(".active").attr("value"));
            $(this).find(".btn-default").html("&nbsp;");
            toggle_best($(this).attr("mode"), $(this).attr("column"),$(this).find(".active").attr("value"));
        });



        $(element + " #species").prop("id", mode_without_space + "_species");
        $(element + " #best_species").prop("id", mode_without_space + "_best_species");
        $(element + " #best_res").prop("id", mode_without_space + "_best_res");

        $(element + " .modal-title").css("display", "inline");
        $(element + " .modal-title").hide();

        if (mode!="Single structure"){
          $(element + " .modal-header").append("<span><button type='button' onclick='check_all(this,1);' class='btn btn-xs btn-primary reset-selection'>Select all displayed</button></span>");
          $(element + " .modal-header").append(" | <span><input type=text class='pastePDBs' placeholder='Paste PDBs with comma- or space-separated'><button type='button' onclick='pastePDBs();' class='btn btn-xs btn-primary reset-selection'>Import</button></span>");
        } else {
          $(element + " .modal-header").append(" <span><input type=text class='pastePDBs' placeholder='Paste PDB'><button type='button' onclick='pastePDBs();' class='btn btn-xs btn-primary reset-selection'>Import</button></span>");
        }

        $(element + " .modal-header .pastePDBs").keypress(function(event) {
            var keycode = (event.keyCode ? event.keyCode : event.which);
            if(keycode === "13"){
               pastePDBs();
            }
        });


        $(element + " .modal-header").append(" | <span><button type='button' onclick='exportPDBs();' class='btn btn-xs btn-primary export_pdbs'>Export</button></span>");
        // $(element + ' .modal-header').append(' | <span><button type="button" onclick="toggle_best(\''+mode+'\',7);" class="btn btn-xs btn-primary">Best</button></span>');
        /*if (window.location.href.endsWith("contactnetwork/clustering") || window.location.href.endsWith("contactnetwork/clustering#"))
          $(element + ' .modal-header').append(' | <span>Structure shortest distance to all other structures of the same receptor and same state: <button type="button" onclick="check_all_distance_representatives();" class="btn btn-xs btn-primary">Distance Representative</button></span>');
        else {
          // a$(element + ' .tableview').before(' | <span>Structure with highest % identity to GPCR’s contact consensus: <button type="button" onclick="check_all_representatives();" class="btn btn-xs btn-primary">Contact Representative</button></span>');
          // $(element + ' .tableview').before(' | <span>Structure sharing either highest/lowest diff between fraction of active/inactive class consensus contacts, or for intermediate the one closes to a 0 diff: <button type="button" onclick="check_all_class_representatives();" class="btn btn-xs btn-primary">New Representative</button></span>');
        }*/
        // $(element + ' .modal-header').append(' | <div class="externalfilters" style="display: inline-block;"><span id="'+mode_without_space+'_external_filter_container_0"></span></div>');
        // $(element + ' .tableview').before('<div class="externalfilters" style="display: inline-block;"><span id="'+mode_without_space+'_external_filter_container_1"></span></div>');

        console.time('DataTable');
        oTable[mode] = $(element + " .tableview table").DataTable({
            "fnInfoCallback": function (oSettings, iStart, iEnd, iMax, iTotal, sPre) {
                filtered = iMax - iEnd;
                filtered_text = filtered ? " (" + filtered + " structures filtered out)" : "";
                var cols = []
                var table = $(element + " .dataTables_scrollBody .structure_selection");
                cols_of_interest = [1, 11];
                for (let [i, row] of [...table.find("tbody")[0].rows].entries()) {
                    for (let [j, cell] of [...row.cells].entries()) {
                        if (cols_of_interest.includes(j)) {
                            cols[j] = cols[j] || [];
                            cols[j].push(cell.innerText)
                        }
                    }
                }
                distinctReceptors = [...new Set(cols[1])];
                distinctReceptorState = [...new Set(cols[1].map((val, i) => [cols[11]].reduce((a, arr) => [...a, arr[i]], [val])))];
                distinctReceptorState = [...new Set(distinctReceptorState.map(x => x[0] + "_" + x[1]))]

                return "Showing " + iEnd + " structures for "+distinctReceptors.length+" receptors and "+distinctReceptorState.length+" distinct receptor-state pairs"+filtered_text;
              },
            "scrollX": true,
            // "paging": true,
            // "autoWidth": true,

            scrollY: "65vh",
            // scrollCollapse: true,
            paging: false,
            "bSortCellsTop": true,
            columnDefs: [{
                targets: "no-sort",
                orderable: false
            },
                {"targets": [ 7, 12 ],
                "visible": false}],
            "aaSorting": [],
            "columns": [
                {
                  "orderDataType": "dom-checkbox"
                },
                null,
                null,
                null,
                null,
                null,
                null,
                null,
                null,
                null,
                null,
                null,
                null,
                null,
                null,
                null,
                null,
                null,
                null,
                null,
                null,
                null,
                null,
                null,
            ],
        });
        console.timeEnd("DataTable");
        console.time("yadcf");
        yadcf.init(oTable[mode],
            [{
                    column_number: 1,
                    filter_type: "multi_select",
                    select_type: "select2",
                    filter_default_label: "UniProt",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 2,
                    filter_type: "multi_select",
                    select_type: "select2",
                    column_data_type: "html",
                    html_data_type: "text",
                    filter_default_label: "GtoPdb",
                    filter_match_mode: "exact",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 3,
                    filter_type: "multi_select",
                    select_type: "select2",
                    html_data_type: "text",
                    select_type_options: {
                        width: "150px"
                    },
                    filter_default_label: "Rec Family",
                    filter_match_mode: "exact",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 4,
                    filter_type: "multi_select",
                    select_type: "select2",
                    filter_default_label: "Cl",
                    select_type_options: {
                        width: "30px"
                    },
                    filter_reset_button_text: false,
                },
                {
                    column_number: 5,
                    filter_type: "range_number",
                    select_type_options: {
                        width: "70px"
                    },
                    filter_default_label: ["From","to"],
                    filter_reset_button_text: false,
                },
                {
                    column_number: 6,
                    filter_type: "multi_select",
                    filter_container_id: mode_without_space + "_species",
                    select_type: "select2",
                    column_data_type: "html",
                    filter_default_label: "Species",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 7,
                    filter_type: "select",
                    select_type: "select2",
                    select_type_options: {
                        width: "90px",
                        minimumResultsForSearch: -1 // remove search box
                    },
                    filter_default_label: "All",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 8,
                    filter_type: "range_number",
                    select_type_options: {
                        width: "70px"
                    },
                    filter_default_label: ["From","to"],
                    filter_reset_button_text: false,
                },
                {
                    column_number: 9,
                    filter_type: "multi_select",
                    select_type: "select2",
                    filter_default_label: "Type",
                    select_type_options: {
                        width: "70px"
                    },
                    filter_reset_button_text: false,
                },
                {
                    column_number: 10,
                    filter_type: "multi_select",
                    select_type: "select2",
                    select_type_options: {
                        width: "70px"
                    },
                    filter_default_label: "PDB",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 11,
                    filter_type: "range_number",
                    select_type_options: {
                        width: "70px"
                    },
                    filter_default_label: ["Res (Å)",""],
                    filter_reset_button_text: false,
                },
                // {
                //     column_number: 12,
                //     filter_type: "multi_select",
                //     select_type: "select2",
                //     filter_default_label: "Best res.",
                //     select_type_options: {
                //         width: "70px"
                //     },
                //     filter_reset_button_text: false,
                // },
                {
                    column_number: 12,
                    // filter_container_id: mode_without_space + "_best_res",
                    filter_type: "select",
                    select_type: "select2",
                    select_type_options: {
                        width: "90px",
                        minimumResultsForSearch: -1 // remove search box
                    },
                    filter_default_label: "All",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 13,
                    filter_type: "multi_select",
                    select_type: "select2",
                    filter_default_label: "State",
                    select_type_options: {
                        width: "70px"
                    },
                    filter_match_mode: "exact",
                    filter_reset_button_text: false,

                },
                // {
                //     column_number : 10,
                //     filter_type: "multi_select",
                //     select_type: "select2",
                //     filter_default_label: "",
                //     filter_match_mode : "exact",
                //     filter_reset_button_text: false,

                // },
                // {
                //     column_number : 11,
                //     filter_type: "multi_select",
                //     select_type: "select2",
                //     filter_default_label: "",
                //     filter_match_mode : "exact",
                //     filter_reset_button_text: false,

                // },
                {
                    column_number: 14,
                    filter_type: "range_number",
                    select_type_options: {
                        width: "70px"
                    },
                    filter_default_label: ["From","to"],
                    // filter_default_label: "Gprot-bound likeness",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 15,
                    filter_type: "range_number",
                    select_type_options: {
                        width: "70px"
                    },
                    filter_default_label: ["From","to"],
                    filter_reset_button_text: false,

                },
                /*{
                    column_number: 13,
                    filter_type: "multi_select",
                    select_type: "select2",
                    filter_default_label: "7TM Open IC (Å)",
                    filter_reset_button_text: false,
                },*/
                {
                    column_number: 16,
                    filter_type: "multi_select",
                    select_type: "select2",
                    filter_default_label: "Family",
                    select_type_options: {
                        width: "70px"
                    },
                    filter_reset_button_text: false,
                },
                {
                    column_number: 17,
                    filter_type: "multi_select",
                    select_type: "select2",
                    filter_default_label: "Subtype",
                    select_type_options: {
                        width: "70px"
                    },
                    filter_reset_button_text: false,
                },
                {
                    column_number: 19,
                    filter_type: "range_number",
                    select_type_options: {
                        width: "50px"
                    },
                    column_data_type: "html",
                    filter_default_label: ["From","to"],
                    filter_reset_button_text: false,
                },
                {
                    column_number: 20,
                    filter_type: "multi_select",
                    select_type: "select2",
                    filter_default_label: "Fusion",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 21,
                    filter_type: "multi_select",
                    select_type: "select2",
                    filter_default_label: "Antibody",
                    column_data_type: "html",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 22,
                    filter_type: "multi_select",
                    select_type: "select2",
                    filter_default_label: "Ligand",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 23,
                    filter_type: "multi_select",
                    select_type: "select2",
                    filter_default_label: "Modality",
                    filter_match_mode: "exact",
                    filter_reset_button_text: false,
                },
                // {
                //     column_number: 24,
                //     filter_type: "multi_select",
                //     select_type: "select2",
                //     filter_default_label: "Modality",
                //     filter_reset_button_text: false,
                // },
                // {
                //     column_number: 25,
                //     filter_container_id: mode_without_space+"_external_filter_container_0",
                //     html_data_type: "text",
                //     select_type: "select2",
                //     // filter_type: "multi_select",
                //     filter_default_label: "All species and structures",
                //     filter_reset_button_text: false,
                //     text_data_delimiter: ",",
                //     select_type_options: {
                //         width: "300px",
                //         minimumResultsForSearch: -1 // remove search box
                //     },
                // },
                // {
                //     column_number: 23,
                //     filter_container_id: mode_without_space+"_external_filter_container_1",
                //     html_data_type: "text",
                //     select_type: "select2",
                //     // filter_type: "multi_select",
                //     filter_default_label: "All Structures",
                //     filter_reset_button_text: false,
                //     select_type_options: {
                //         width: "250px",
                //         minimumResultsForSearch: -1 // remove search box
                //     },
                // },
            ], {
                cumulative_filtering: false,
                // filters_tr_index: 1
            }
        );

        console.timeEnd("yadcf");
        oTable[mode].columns.adjust();

        // console.time("yadcf_reset");
        // yadcf.exResetAllFilters(oTable[mode]);
        // console.timeEnd("yadcf_reset");
        // yadcf.exFilterColumn(oTable[mode], [
        //     [21, "*Only show mammalian structures and those from human or closest species"],
        //   ]);

        var loading_anim = $(element + " .loading_overlay");
        var delayDraw = true;

        oTable[mode].on('preDraw.dt', function() {
            console.time('preDrawing');
            console.timeEnd('preDrawing')
            if (delayDraw) {
                console.log('ON ---------------------------------');
                console.time('showing');
                loading_anim.show();
                console.timeEnd('showing');
                console.log(loading_anim.css('display'));

                setTimeout(function() {
                delayDraw = false;
                oTable[mode].draw();
                }, 100);

                return false;
            }
        });

        oTable[mode].on("draw.dt", function(e, oSettings) {
            console.time("create_overlay");
            create_overlay(element + " .structure_selection");
            console.timeEnd("create_overlay");
            loading_anim.hide();
            delayDraw = true;
            console.time("update_text_in_modal");
            update_text_in_modal();
            console.timeEnd("update_text_in_modal");
        });




        // $(element + ' .dataTables_scrollBody').append('<div class="structure_overlay"><table id="overlay_table" class="overlay_table row-border text-center compact dataTable no-footer text-nowrap"><tbody></tbody></table></div>');
        $(element + ' .dataTables_scroll').append('<div class="structure_overlay"><table id="overlay_table" class="overlay_table row-border text-center compact dataTable no-footer text-nowrap"><tbody></tbody></table></div>');
        $(element + " .structure_overlay").hide();

        $(element + ' .dataTables_scrollBody').before("<div class='top_scroll'><div>&nbsp;</div></div>");

        var dataTables_scrollBody_height = $(element + " .dataTables_scrollBody")[0].offsetHeight;
        var bodyRect = $(element + " .dataTables_scroll")[0].getBoundingClientRect(),
        elemRect = $(element + " .dataTables_scrollBody")[0].getBoundingClientRect(),
        offset = elemRect.top - bodyRect.top;

        $(".top_scroll").css({
            "width": "100%",
            "overflow-x": "scroll",
            "overflow-y": "auto",
        });
        $(".top_scroll div").css({
            // "background-color": "red",
            "font-size": "1px",
            "line-height": "1px",
        });

        $(".structure_overlay").css({
            // "top": "0px",
            "top": offset+"px",
            "position": "absolute",
            "background": "#f8f8f8",
            "-webkit-box-shadow": "5px 0 2px -2px #888",
            "box-shadow": "5px 0 2px -2px #888",
            "height": (dataTables_scrollBody_height-17) + "px",
            "overflow-y": "scroll",
            "scrollbar-width": "none",
        });


        $(".structure_overlay tbody tr").css({
            "background-color": "#f8f8f8",
        });

        create_overlay(element + " .structure_selection");
        // track_scrolling(element);
        new_track_scrolling(element);
        console.log("done overlays");


        $(function() {
            var tableContainer = $(".dataTables_scrollBody");
            var table = $(".dataTables_scrollBody table");
            var fakeContainer = $(".top_scroll");
            var fakeDiv = $(".top_scroll div");

            var tableWidth = table.width();
            fakeDiv.width(tableWidth);

            fakeContainer.scroll(function() {
                tableContainer.scrollLeft(fakeContainer.scrollLeft());
            });
            tableContainer.scroll(function() {
                fakeContainer.scrollLeft(tableContainer.scrollLeft());
            });
        })
        console.log("done scrolling");

        $(".dataTables_scrollBody tr").css("cursor", "pointer");
        $(".dataTables_scrollBody tr").click(function(event) {
            if (event.target.type !== "checkbox") {
                $(":checkbox", this).trigger("click");
                if ($(":checkbox", this).length === 0) {
                    var pdb_id = $(this).attr("id").split("_")[1];
                    var checkbox = $(this).closest(".dataTables_scrollBody").find("#" + pdb_id);
                    checkbox.trigger("click");
                }
            }
        });
        console.log("done click event");

    };
    $(element + " .loading_overlay").hide();
    prepopulatePDBs();
}

function track_scrolling(element) {
    var left = 0;
    var old_left = 0;
    toggle_enabled = true;
    var d = new Date();
    last_change = d.getTime();
    $(element + " .dataTables_scrollBody").scroll(function () {

        var d = new Date();
        var now = d.getTime();

        if (now>(last_change+1000)) {

            // If user scrolls and it"s >100px from left, then attach fixed columns overlay
            left = $(element + " .dataTables_scrollBody").scrollLeft();
            // if (left != old_left) $(".structure_overlay").hide();
            // old_left = left;


            console.log(now, last_change);

            if (left > 100 && toggle_enabled && left != old_left) {
                console.log("SCROLL ", left, now);
                $(".structure_overlay").css({
                    left: left + "px"
                });
                if ($(".structure_overlay").is(":hidden")){
                    // link check boxes
                    $(".structure_overlay").show();
                }
            } else if (left < 100) {
                if ($(".structure_overlay").is(":visible")){
                // link check boxes
                $(".structure_overlay").hide();
                }
            }
            old_left = left;
            last_change = now;
        }
    });
}

function new_track_scrolling(element) {
    var left = 0;
    var old_left = 0;
    toggle_enabled = true;
    old_height = 0;
    old_top = 0;
    var d = new Date();
    var scrollTop;
    last_change = d.getTime();
    scroll_visible = false;

    // Init the size to begin with.. It is updated once in a while incase of window resizing.
    var dataTables_scrollBody_height = $(element + " .dataTables_scrollBody")[0].offsetHeight;
    var bodyRect = $(element + " .dataTables_scroll")[0].getBoundingClientRect(),
        elemRect = $(element + " .dataTables_scrollBody")[0].getBoundingClientRect(),
        offset = elemRect.top - bodyRect.top;
    $(".structure_overlay").css({
        height: (dataTables_scrollBody_height - 17) + "px",
        top: offset + "px"
    });

    $(element + " .dataTables_scrollBody").scroll(function () {
        var d = new Date();
        now = d.getTime();

        // This isnt dom heavy..
        scrollTop = $(element + " .dataTables_scrollBody").scrollTop();
        if (scrollTop != old_top) {
            $(".structure_overlay").scrollTop(scrollTop);
            old_top = scrollTop
        }

        left = $(element + " .dataTables_scrollBody").scrollLeft();
        if ( left<100 && scroll_visible){
            $(".structure_overlay").hide();
            scroll_visible = false;
        } else if ( left>100 && !scroll_visible){
            $(".structure_overlay").show();
            scroll_visible = true;
        }

        if (now > (last_change + 1000)) {

            // If user scrolls and it"s >100px from left, then attach fixed columns overlay
            if (left > 100 && toggle_enabled) {

                var dataTables_scrollBody_height = $(element + " .dataTables_scrollBody")[0].offsetHeight;
                if (dataTables_scrollBody_height != old_height) {
                    var bodyRect = $(element + " .dataTables_scroll")[0].getBoundingClientRect(),
                        elemRect = $(element + " .dataTables_scrollBody")[0].getBoundingClientRect(),
                        offset = elemRect.top - bodyRect.top;
                    $(".structure_overlay").css({
                        height: (dataTables_scrollBody_height - 17) + "px",
                        top: offset + "px"
                    });
                    old_height = dataTables_scrollBody_height;
                }
                if (!scroll_visible) {
                    $(".structure_overlay").show();
                    scroll_visible = true;
                }
            }
            last_change = now;
        }
    });
    $(".structure_overlay").scroll(function () {
        scrollTop = $(".structure_overlay").scrollTop();
        if (scrollTop != old_top) {
            $(element + " .dataTables_scrollBody").scrollTop(scrollTop);
            old_top = scrollTop
        }
    });
}
