/*global yadcf, rankedRangeFilter, createRank, createYADCFfilters, supportFilter, make_rank_col_filters, lastRangeRankFilter, initFixedColumnsOverlay, copyListToClipboard*/
/*eslint complexity: ["error", 20]*/
/*eslint wrap-iife: ["error", "outside"]*/
/*eslint quotes: ["error", "double", { "avoidEscape": true }]*/

let oTable1 = [];

function rankedRangeFiltert1(filterVal, columnVal, rowValues, stateVal) {
  return rankedRangeFilter(filterVal, columnVal, rowValues, stateVal, "3");
}

$(document).ready(function() {
  // Activate tooltips and popovers from Bootstrap
  $("[data-toggle='tooltip']").tooltip();
  $("[data-toggle='popover']").popover();

  // =============================================================================
  // Families Table
  // =============================================================================
  console.time("table1load");
  $("#arrestintable").show();
  oTable1 = $("#arrestintable").DataTable({
    deferRender: true,
    scrollY: '55vh',
    scrollX: true,
    scrollCollapse: true,
    scroller: true,
    paging: false,
    bSortCellsTop: false, //prevent sort arrows going on bottom row
    aaSorting: [],
    order: [
      [2, "asc"],
      [4, "asc"]
    ],
    autoWidth: true,
    bInfo: true,
  });

  let column_filters = [];
  // Selector column
  column_filters = column_filters.concat(createYADCFfilters(0, 1, "none"));
  // Receptor section
  column_filters = column_filters.concat(createYADCFfilters(1, 1, "multi_select", "select2", "", false, null, "html", "80px"));
  column_filters = column_filters.concat(createYADCFfilters(2, 1, "multi_select", "select2", "", false, null, null, "40px"));
  column_filters = column_filters.concat(createYADCFfilters(3, 1, "multi_select", "select2", "", false, null, null, "200px"));
  column_filters = column_filters.concat(createYADCFfilters(4, 1, "multi_select", "select2", "", false, "exact", "html", "80px"));
  column_filters = column_filters.concat(createYADCFfilters(5, 1, "multi_select", "select2", "", false, "exact", "html", "80px"));
  // Ligands section
  column_filters = column_filters.concat(createYADCFfilters(6, 1, "multi_select", "select2", "", false, "exact", "html", "200px"));
  column_filters = column_filters.concat(createYADCFfilters(7, 1, "multi_select", "select2", "", false, null, null, "50px"));
  // Coupling data filters
  column_filters = column_filters.concat(createYADCFfilters(8, 9, "range_number", null, ["Min", "Max"], false, null, null, "40px", null, "-"));
  // column_filters = column_filters.concat(make_rank_col_filters(8, 19, "hide_ranksub", rankedRangeFiltert1));
  // Hidden GPCRdb support type column calls customized function
  // column_filters = column_filters.concat([
  //   {
  //     column_number: 21,
  //     filter_type: "custom_func",
  //     custom_func: supportFilter,
  //     filter_container_id: "hide_filter3",
  //   }]);

  yadcf.init(oTable1, column_filters, {
    cumulative_filtering: false
  });

  // Initialize ranked Range Filtering options
  // $(".ranked_range_min3, .ranked_range_max3, .ranked_range_rank3").on("click", function(event) {
  //   event.stopPropagation();
  // });

  // $(".ranked_range_min3, .ranked_range_max3, .ranked_range_rank3").on("input", function(event) {
  //   // Get column #
  //   let column_nr = event.target.id.split("_")[3];

  //   // Store current type of filtering globally
  //   lastRangeRankFilter = event.target.id.split("_")[2];

  //   // SELECT one option in real YADCF filter to trick YADCF into calling the filter function
  //   let adjust_node1 = $("#yadcf-filter--arrestintable-" + column_nr + " option").filter(":nth-child(2)").first();
  //   adjust_node1.prop("selected", true);

  //   // Invoke filtering
  //   $("#yadcf-filter--arrestintable-" + column_nr).change();

  //   // Clean filter type
  //   lastRangeRankFilter = "";
  // });

  // Reset filters + recalculates table layout
  // yadcf.exResetAllFilters(oTable1);

  $("#arrestintable" + " > tbody > tr").click(function(event) {
    if (event.target.type !== "checkbox") {
      $(":checkbox", this).trigger("click");
      $(this).eq(0).toggleClass("alt_selected");
      $(this).find("td").toggleClass("highlight");
    }
    $(this).eq(0).toggleClass("alt_selected");
    $(this).find("td").toggleClass("highlight");
  });
  console.timeEnd("table1load");

  $(".select-all").click(function() {
    $(":checkbox", this).trigger("click");
    if ($(this).prop("checked") === true) {
      $(".alt").prop("checked", true);
      $(".alt").parent().parent().addClass("alt_selected");
      $(".alt").parent().parent().find("td").addClass("highlight");
    }
    if ($(this).prop("checked") === false) {
      $(".alt").prop("checked", false);
      $(".alt").parent().parent().removeClass("alt_selected");
      $(".alt").parent().parent().find("td").removeClass("highlight");
    }
  });

  // Put top scroller
  // https://stackoverflow.com/questions/35147038/how-to-place-the-datatables-horizontal-scrollbar-on-top-of-the-table
  //    console.time("scroll to top");
  // $(".dataTables_scrollHead").css({
  //   "overflow-x": "scroll"
  // }).on("scroll", function(e) {
  //   var scrollBody = $(this).parent().find(".dataTables_scrollBody").get(0);
  //   scrollBody.scrollLeft = this.scrollLeft;
  //   $(scrollBody).trigger("scroll");
  // });

  // Enable columns overlay
  initFixedColumnsOverlay("arrestintable");

  // Hide column buttons
  $(".hide_columns3").click(function(evt) {
    var columns = $(this).attr("columns").split(",");
    oTable1.columns(columns).visible(false, false);
    oTable1.draw();
    oTable1.columns.adjust();
  });

  // Enable copy to clipboard option
  $("#uniprot_copy1").click(function() {
    let selected = [];
    $(".alt_selected > .uniprot1 > a").each(function() {
      selected.push($(this).text());
    });
    copyListToClipboard(selected);
  });

  //Uncheck every row when using back button on browser
  $(".alt_selected").prop("checked", false);
  $(".alt").prop("checked", false);
  $(".select-all").prop("checked", false);

  oTable1.columns.adjust();

});
