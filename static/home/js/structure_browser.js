/*eslint complexity: ["error", 8]*/
/*global ClearSelection, CheckSelection, AddToSelection, superposition, copyToClipboard, showAlert*/

function structurebrowser() {

    var oTable2 = $("#structures_scrollable").DataTable({
        scrollY:        "65vh",
        scrollX:        true,
        scrollCollapse: true,
        scroller:       true,
        paging:         true,
        pageLength:     100, 
        aaSorting:      [],
        autoWidth:      false,
        order:          [[29, "desc"],[1,"asc"]],
        columnDefs:      [{ "targets": "no-sort", "orderable": false }],
        columns:        [
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
                          {"width": "20%"},
                          null,
                          null,
                          null,
                          null,
                          null,
                          null,
                          null,
                          null,
                          null,
                          null
                      ],
        bInfo:        true,
    });

    var prev_ids = Array();
    var current_align_ids = Array();

    //Uncheck every row when using back button on browser
    $(".alt_selected").prop("checked",false);
    $(".alt").prop("checked",false);
    $(".select-all").prop("checked",false);

    const validPrefixes = ["#keepselection"];

    if (validPrefixes.some(prefix => window.location.hash.startsWith(prefix))) {}
    else {
      ClearSelection('targets');
      ClearSelection('reference');
    }

    $("#loading_div").hide();

    yadcf.init(oTable2,
        [
            {
                column_number : 1,
                filter_type: "multi_select",
                select_type: "select2",
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "UniProt",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "60px",
                }
            },
            {
                column_number : 2,
                filter_type: "multi_select",
                select_type: "select2",
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "IUPHAR",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "60px",
                }
            },
            {
                column_number : 3,
                filter_type: "multi_select",
                select_type: "select2",
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "Select",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "80px",
                }
            },
            {
                column_number: 4,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "35px",
                }
            },
            {
                column_number: 5,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "Species",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "60px",
                }
            },
            {
                column_number : 6,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "Method",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "60px",
                }
            },
            {
                column_number : 7,
                filter_type: "multi_select",
                select_type: "select2",
                column_data_type: "html",
                filter_default_label: "Select",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "50px",
                }
            },
            {
                column_number : 8,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "Select",
                column_data_type: "html",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "50px",
                }
            },
            {
                column_number : 9,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 10,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "30px",
                }
            },
            {
                column_number : 11,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "State",
                filter_reset_button_text: false,
                filter_match_mode : "exact",
                select_type_options: {
                    width: "70px",
                }
            },
            {
                column_number : 12,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 13,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: "50px",
                }
            },
            {
                column_number : 14,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "Family",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "50px",
                }
            },
            {
                column_number : 15,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "Subtype",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "55px",
                }
            },
            {
                column_number : 16,
                filter_type: "text",
                select_type: "select2",
                filter_default_label: "Note",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "80px",
                }
            },
            {
                column_number : 17,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: "50px",
                }
            },
            {
                column_number : 18,
                filter_type: "text",
                select_type: "select2",
                filter_default_label: "Fusion",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "100px",
                }
            },
            {
                column_number : 19,
                filter_type: "text",
                select_type: "select2",
                filter_default_label: "Antibodies",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "100px",
                }
            },
            {
                column_number : 20,
                filter_type: "text",
                select_type: "select2",
                html_data_type: "text",
                filter_default_label: "Ligand name",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "100px",
                }
            },
            {
                column_number : 21,
                filter_type: "multi_select",
                select_type: "select2",
                html_data_type: "text",
                filter_default_label: "Ligand type",
                filter_reset_button_text: false,
                text_data_delimiter: "<br>",
                select_type_options: {
                    width: "100px",
                }
            },
            {
                column_number : 22,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "Modality",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                text_data_delimiter: "<br>",
                select_type_options: {
                    width: "100px",
                }
            },
            {
                column_number : 23,
                filter_type: "text",
                select_type: "select2",
                filter_default_label: "Ligand name",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "100px",
                }
            },
            {
                column_number : 24,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "Ligand type",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "80px",
                }
            },
            {
                column_number : 25,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "40px",
                }
            },
            {
                column_number : 26,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "40px",
                }
            },
            {
                column_number : 27,
                filter_type: "text",
                select_type: "select2",
                filter_default_label: "Authors",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "100px",
                }
            },
            {
                column_number : 28,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "Reference",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "140px",
                }
            },
            {
                column_number : 29,
                filter_type: "range_date",
                filter_reset_button_text: false,
                date_format: "yyyy-mm-dd",
                select_type: "select2",
                filter_default_label: ["Min", "Max"],
                // filter_reset_button_text: false,
            },
            {
                column_number : 30,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "Select",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "80px",
                }
            }
        ],
        {
            cumulative_filtering: false
        }
    );

    // NOTE: getting the columns correct with DT/YADCF
    //yadcf.exResetAllFilters(oTable2); // Slowest option
    //oTable2.draw();                   // Bit faster but still slow
    oTable2.columns.adjust()            // Quick and clean

    $("#structures_scrollable > tbody > tr").click(function(event) {
        if (event.target.type !== "checkbox") {
            $(":checkbox", this).trigger("click");
            $(this).eq(0).toggleClass("alt_selected");
            $(this).find("td").toggleClass("highlight");
        }
        $(this).eq(0).toggleClass("alt_selected");
        $(this).find("td").toggleClass("highlight");
    });

    $(".select-all").click(function() {
        $(":checkbox", this).trigger("click");
        if ($(this).prop("checked")===true) {
            $(".alt").prop("checked", true);
            $(".alt").parent().parent().addClass("alt_selected");
            $(".alt").parent().parent().find("td").addClass("highlight");
        }
        if ($(this).prop("checked")===false) {
            $(".alt").prop("checked", false);
            $(".alt").parent().parent().removeClass("alt_selected");
            $(".alt").parent().parent().find("td").removeClass("highlight");
        }
    });

    $(".hide_columns").click(function(evt) {
    var columns = $(this).attr("columns").split(",");
    columns.forEach(function(column) {
        var column = oTable2.column( column );
        try {
            column.visible( false, false );
        }
        catch(err) {
            column.visible( false, false );
        }
        });
        oTable2.draw();
    } );

    var toggle_enabled = true;
    $("#toggle_fixed_btn").click(function() {
        if (toggle_enabled) {
            toggle_enabled = false;
            $("#overlay").hide();
            $("#toggle_fixed_btn").attr("value","Enable fixed columns");
            $("#toggle_fixed_btn").addClass("clicked_button");
        } else {
            toggle_enabled = true;
            $(".dataTables_scrollBody").scroll();
            $("#toggle_fixed_btn").attr("value","Disable fixed columns");
            $("#toggle_fixed_btn").removeClass("clicked_button");
        }
    });

    $("#toggle_columns_btn").click(function() {
        var columns = Array.from(new Array(24), (x,i) => i + 6);
        columns.forEach(function(column) {
            var column = oTable2.column( column );
            try {
                column.visible( true, false );
            }
            catch(err) {
                column.visible( true, false );
            }
        });
        oTable2.draw();
    });

    $("#representative_btn").click(function () {
        $(this).toggleClass("toggled");
        if ($(this).hasClass("toggled")) {
            $("#representative_btn").addClass("clicked_button");
            $("#representative_btn").attr("value","All structures");
        }
        else {
            $("#representative_btn").removeClass("clicked_button");
            $("#representative_btn").attr("value","Representative structures (state & receptor)");
        }
        oTable2.draw();
    });

    $.fn.dataTable.ext.search.push(
        function (settings, data, dataIndex) {
            if ($("#representative_btn").hasClass("toggled")) {
                if ($(oTable2.row(dataIndex).node()).hasClass("repr-st") && $(oTable2.row(dataIndex).node()).hasClass("repr-st")) {
                    return true;
                }
                return false;
            }
            else {
                return true;
            }
    });

    $("#align_btn").click(function () {
        var checked_data = oTable2.rows(".alt_selected").data();
        ClearSelection("targets");
        if (checked_data.length === 0) {
            showAlert("No entries selected for sequence alignment", "danger");
            return 0;
        }
        for (var i = 0; i < checked_data.length; i++) {
            var div = document.createElement("div");
            div.innerHTML = checked_data[i][7];
            if (typeof div.innerText !== "undefined") {
                AddToSelection("targets", "structure", div.innerText.replace(/\s+/g, ""));
            } else {
                AddToSelection("targets", "structure", div.textContent.replace(/\s+/g, ""));
            }
        }
        window.location.href = "/structure/selection_convert";
    });

    $("#superpose_btn").click(function() {
        superposition(oTable2, [7,1,2,3,4,5,11,29], "structure_browser", "gpcr", 7);
    });

    $('#superpose_template_btn').click(function () {
        direct_superposition(oTable2, "gpcr", 7, "structure");
    });

    $("#download_btn").click(function () {
        ClearSelection("targets");
        var checked_data = oTable2.rows(".alt_selected").data();
        if (checked_data.length === 0) {
            showAlert("No structures selected for download", "danger");
            return 0;
        }
        else if (checked_data.length > 100) {
            showAlert("Maximum number of selected entries is 100", "warning");
            return 0;
        }
        var selected_ids = [];
        for (var i = 0; i < checked_data.length; i++) {
            var div = document.createElement("div");
            div.innerHTML = checked_data[i][7];
            if (typeof div.innerText !== "undefined") {
                selected_ids.push(div.innerText.replace(/\s+/g, ""));
            } else {
                selected_ids.push(div.textContent.replace(/\s+/g, ""));
            }
        }
        AddToSelection("targets", "structure_many", selected_ids.join(","));

        window.location.href = "/structure/pdb_download";
    });

    $(".uniprot-export").data("powertipjq", $([
        "<p>Export UniProt IDs</p>"
        ].join("\n")));
    $(".pdb-export").data("powertipjq", $([
        "<p>Export PDB IDs</p>"
        ].join("\n")));
    $(".glyphicon-export").powerTip({
        placement: "n",
        smartPlacement: true
    });
    $("#uniprot_copy").click(function () {
        copyToClipboard($(".alt_selected > .uniprot > a"), "\n", "UniProt IDs", $(".uniprot-export"));
    });
    $("#pdb_copy").click(function () {
        copyToClipboard($(".alt_selected > .pdb > a"), "\n", "PDB IDs", $(".pdb-export"));
    });

    $("#reset_filters_btn").click(function () {
        window.location.href = "/structure/";
    });

    $(".dataTables_scrollBody").append("<div id=overlay><table id='overlay_table' class='row-border text-center compact dataTable no-footer text-nowrap'><tbody></tbody></table></div>");

    function create_overlay() {
        // This function fires upon filtering, to update what rows to show as an overlay
        $("#overlay_table tbody tr").remove();
        var $target = $("#overlay_table tbody");
        $("#structures_scrollable tbody tr").each(function() {
            var $tds = $(this).children(),
                $row = $("<tr></tr>");
            // $row.append($tds.eq(0).clone()).append($tds.eq(1).clone()).appendTo($target);
            $row.append($tds.eq(1).clone()).append($tds.eq(2).clone()).appendTo($target);
            $row.height($(this).height());
        });
        $("#overlay_table .border-right").removeClass("border-right");
    }

    // Function that detects filtering events
    $("#structures_scrollable").on( "draw.dt", function (e,oSettings) {
        create_overlay();
    });

    create_overlay();
    $("#overlay").hide();

    var left = 0;
    var old_left = 0;
    $(".dataTables_scrollBody").scroll(function(){
        // If user scrolls and it"s >100px from left, then attach fixed columns overlay
        left = $(".dataTables_scrollBody").scrollLeft();
        if (left !== old_left) {
          $("#overlay").hide();
        }
        old_left = left;

        if (left>100 && toggle_enabled) {
            $("#overlay").css({ left: left+"px" });
            if ($("#overlay").is(":hidden")) {
              $("#overlay").show();
            }
        }
    });
}
