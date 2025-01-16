function showColumn(colClass, event) {
    // Hide all columns
    $('.icl2, .icl3, .c-term, .icl2_icl3, .icl2_c-term, .icl3_c-term, .icl2_icl3_c-term').hide();
    // $('.icl2, .icl3, .c-term').hide();
    console.log('##### showColumn ######')

    // Show the selected column
    $('.' + colClass).show();

}

// Truncate text content in elements with class 'truncate-js'

$(document).ready(function() {

    $('.truncate').each(function() {
        let fullText = $(this).text();
        // console.log('Original Text:', fullText);
    
        if (fullText.length > 3) {
            $(this).text(fullText.substring(0, 4) + '...');
            $(this).attr('title', fullText);
        }
    });

    let oTable = $('#table').DataTable({
        deferRender: true,
        scrollY: '65vh',
        scrollX: true,
        scrollCollapse: true,
        scroller: true,
        paging: false,
        bSortCellsTop: false, //prevent sort arrows going on bottom row
        aaSorting: [],
        order: [
          [5, "asc"],
          [1, "asc"]
        ],
        autoWidth: true,
        bInfo: true,
    });
    yadcf.init(oTable, [
            {column_number: 1, filter_type: 'multi_select', select_type:'select2'},
            {column_number: 2, filter_type: 'multi_select', select_type:'select2'},
            {column_number: 3, filter_type: 'multi_select', select_type:'select2'},
            {column_number: 4, filter_type: 'multi_select', select_type:'select2'},
            {column_number: 6, filter_type: 'range_number', filter_default_label:['Min', 'Max']}
    ]);


    // First tab initiation 
    showColumn('icl2');
    // Changing span of phosphorylation cell to match the tabs
    $('.main').attr('colspan', '10')
    oTable.draw()

    // Refresh the table with the current tab
    oTable.on('draw.dt', function() {
        console.log('########## Re drawing')
        let activeColumn = $('.nav-item.active').attr('data-column');
        console.log(activeColumn , 'Current collumn')
        showColumn(activeColumn)
    });

    $("#table" + " > tbody > tr").click(function(event) {
        if (event.target.type !== "checkbox") {
          $(":checkbox", this).trigger("click");
          $(this).eq(0).toggleClass("alt_selected");
          $(this).find("td").toggleClass("highlight");
        }
        $(this).eq(0).toggleClass("alt_selected");
        $(this).find("td").toggleClass("highlight");
      });
});

