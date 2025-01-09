////// BEZIER THINGS ///////
bezier_cache = [];
function wherebezier(p04,p14,p24,step,stop,p34,allow_cache) {
    //https://en.wikipedia.org/wiki/B%C3%A9zier_curve
    pos = 0;
    length = 0;
    var p = p04;
    xy = [0,0];

    i = 0;
    while (pos <= 1) {
        if (length>stop) { //stop if it reached the length along the line
            break;
        }
        if (i in bezier_cache && allow_cache) {
            xy = bezier_cache[i][0];
            length = bezier_cache[i][1];
        } else {

            if (p34 === undefined) {
                xy = bezier(p04,p14,p24,pos);
            } else {
                xy = bezier_high(p04,p14,p24,p34,pos);
            }
            length += Math.round(100*Math.sqrt( Math.pow(xy[0]-p[0],2) + Math.pow(xy[1]-p[1],2) ))/100;
        }
        p = xy;
        pos += step;
        bezier_cache[i] =  [xy,length]
        i += 1;
    }
    return [xy,length];
}

function lengthbezier(p03,p13,p23,step,p33) {
    //https://en.wikipedia.org/wiki/B%C3%A9zier_curve

    pos = 0;
    length = 0;
    p = p03;
    while (pos <= 1) {

        if (p33 === undefined) {
            xy = bezier(p03,p13,p23,pos);
        } else {
            xy = bezier_high(p03,p13,p23,p33,pos);
        }

        length += Math.sqrt( Math.pow(xy[0]-p[0],2) + Math.pow(xy[1]-p[1],2) );
        p = xy;
        pos += step;
    }
    return Math.round(length);
}

function bezier(p02,p12,p22,t){
    //https://en.wikipedia.org/wiki/B%C3%A9zier_curve
    v1x = p12[0]-p02[0];
    v1y = p12[1]-p02[1];

    i12 = [p02[0]+(p12[0]-p02[0])*t,p02[1]+(p12[1]-p02[1])*t];
    i22 = [p12[0]+(p22[0]-p12[0])*t,p12[1]+(p22[1]-p12[1])*t];

    return [i12[0]+(i22[0]-i12[0])*t,i12[1]+(i22[1]-i12[1])*t];
}

function bezier_high(p01,p11,p21,p31,t){
    //https://en.wikipedia.org/wiki/B%C3%A9zier_curve
    i1 = bezier(p01,p11,p21,t);
    i2 = bezier(p11,p21,p31,t);

    return [i1[0]+(i2[0]-i1[0])*t,i1[1]+(i2[1]-i1[1])*t];
}

function add_box(x,y,name,color){
    box_text = "";
    name = name.substring(0, 6);
    n = name.length;
    width = n*7+2;
    box_text += "<rect class='"+name+" construct_custom' x="+(x-width/2)+" y="+(y-10)+" rx=5 ry=5 width='"+width+"' height='20' stroke='"+color+"' fill='white' stroke-width='2' style2='fill:red;stroke:black;stroke-width:5;opacity:0.5'/>";
    box_text += "<text  class='"+name+" construct_custom rtext' x="+(x)+" y="+(y+4)+" text-anchor='middle' font-size=12 font-family='helvetica'>"+name+"</text>";
    return box_text;
}
//////////////////////////



presetColors = {'D': ['#E60A0A', '#FDFF7B'],'E': ['#E60A0A', '#FDFF7B'],
                'K': ['#145AFF', '#FDFF7B'],'R': ['#145AFF', '#FDFF7B'],
                'S': ['#A70CC6', '#FDFF7B'],'T': ['#A70CC6', '#FDFF7B'],
                'N': ['#A70CC6', '#FDFF7B'],'Q': ['#A70CC6', '#FDFF7B'],
                'V': ['#E6E600', '#000000'],'L': ['#E6E600', '#000000'],
                'I': ['#E6E600', '#000000'],'A': ['#E6E600', '#000000'],
                'M': ['#E6E600', '#000000'],'F': ['#18FF0B', '#000000'],
                'Y': ['#18FF0B', '#000000'],'W': ['#0BCF00', '#000000'],
                'H': ['#0093DD', '#000000'],'P': ['#CC0099', '#FDFF7B'],
                'C': ['#B2B548', '#000000'],'G': ['#FF00F2', '#000000'],
                '-': ['#FFFFFF', '#000000'],
                '+': ['#FFFFFF', '#000000']
                };

var translateOffset = 0;
function showToolTip(x, y, str,rid,plotid) {
    var tipElement = document.getElementById('tool-tip-'+plotid);

    var originalCircle = $('#'+plotid).find("#"+rid);

    //console.log(originalCircle.attr('extra'));
    //originalCircle.attr('extra');

    var rect = tipElement.childNodes[1];
    var text = tipElement.childNodes[3];

    while (text.lastChild) {
       text.removeChild(text.lastChild);
    }

    // var NS = "http://www.w3.org/2000/svg";


    // //text.textContent =  str;
    //     var text_tspan = document.createElementNS(NS, "tspan");


    //     rect.setAttribute('height', 25);
    //     rect.setAttribute('y', -40);
    //     text_tspan.textContent = String(str);
    //     text.appendChild(text_tspan);

    if (originalCircle.attr('extra')) { //Only display if mutateddata flag is on and there is info

        text.innerHTML =  "<tspan  x=\"0\" y=\"-33\">" + str + "</tspan>";
        //text.innerHTML =  text.innerHTML + "<tspan  x=\"0\" y=\"-20\">" + residueColor.mutatedCalc[rid][0] + " mutations" + " | " + residueColor.mutatedCalc[rid][1] + " maxFold"+ " | " + residueColor.mutatedCalc[rid][2] + " minFold </tspan>";
        text.innerHTML =  text.innerHTML + "<tspan  x=\"0\" y=\"-20\">" + originalCircle.attr('extra') + " </tspan>";
        rect.setAttribute('height', 35);
        rect.setAttribute('y', -50);
    } else {
        text.innerHTML =  "<tspan  x=\"0\" y=\"-23\">" + str + "</tspan>";
        rect.setAttribute('height', 25);
        rect.setAttribute('y', -40);
    }

    var bbox = text.getBBox();
    rect.setAttribute('width', bbox.width + 8);
    rect.setAttribute('x', -bbox.width/2 - 4);

    var transX = (x <= (bbox.width + 8) / 2) ? (bbox.width + 8) / 2 : x;
    tipElement.setAttribute('transform', 'translate(' + transX + ',' + (y + translateOffset) + ')');
    tipElement.setAttribute('visibility', 'visible');
}

function hideToolTip(plotid) {
    var tipElement = document.getElementById('tool-tip-'+plotid);
    tipElement.setAttribute('visibility', 'hidden');
}

function toggleLoop(id, type, skipmaxmin, el) {
  svg = $(el).closest('svg');
    if (type=='long') {
      svg.find(id+".long").hide();
      svg.find(id+".short").show();
    } else {
      svg.find(id+".long").show();
      svg.find(id+".short").hide();
    }
    // $(id+".long").each(function () {
    //     curr = $(this).css("display");
    //     if (curr == 'none') $(this).removeAttr("display");
    //     if (!curr) $(this).css("display", "none");
    // });

    // $(id+".short").each(function () {
    //     curr = $(this).css("display");
    //     if (curr == 'none') $(this).removeAttr("display");
    //     if (!curr) $(this).css("display", "none");
    // });
    // console.log('draw c', (+new Date()) - start);
/*    if (id == ".ICL3" || id == ".ICL2" || id == ".ICL1" || id == ".C-term" ) {
      if (skipmaxmin!=1) redraw_terminus("C");
    }
    // console.log('draw n', (+new Date()) - start);
    if (id == ".ECL3" || id == ".ECL2" || id == ".ECL1" || id == ".N-term" ) {
      if (skipmaxmin!=1) redraw_terminus("N");
    }*/
    // console.log('maxmin', (+new Date()) - start);
    // if (skipmaxmin!=1) maxmin();
    // console.log('done', (+new Date()) - start);
    // if (skipmaxmin!=1) $("#snake").html($("#snake").html());
    // if (skipmaxmin!=1) reload_tooltips();

    // console.log('toggleLoop');
    if (skipmaxmin!=1 && typeof overlay_modifications === 'function') overlay_modifications();
    if (skipmaxmin!=1) maxmin();
}

function redraw_terminus(term) {
  // console.log('redraw term',term);
  if (term=='C') {
    orientation = 1;
    y_max = 0;
    // y_min = 0;
    var x_max = 0;
    $( "circle.ICL1:visible,circle.ICL2:visible,circle.ICL3:visible" ).each(function() {
      this_y = parseInt($(this).attr('cy'));
      if (this_y>y_max) y_max = this_y;
    });
    $( "rect.ICL1:visible,rect.ICL2:visible,rect.ICL3:visible" ).each(function() {
      this_y = parseInt($(this).attr('y'));
      if (this_y>y_max) y_max = this_y;
    });
    y_max = y_max+50;
  } else {
    // N -term
    orientation = -1;
    y_min = 0;
    // $( "circle.ECL1:visible,circle.ECL2:visible,circle.ECL3:visible" ).each(function() {
    //   this_y = parseInt($(this).attr('cy'));
    //   if (this_y<y_min) y_min = this_y;
    // });
    $( "rect.ECL1:visible,rect.ECL2:visible,rect.ECL3:visible" ).each(function() {
      this_y = parseInt($(this).attr('y'));
      if (this_y<y_min) y_min = this_y;
    });
    y_min = y_min-50;
    var x_max = x_svgmax; //from minmax function
  }
    // console.log('maxmin found');

  // var x1 =
  var residues = [];
  // generate list of positions to move
  $( "circle."+term+"-term:visible" ).each(function() {
      id = $(this).attr('id');
      residues.push(parseInt(id));
      // this_y = $(this).attr('y');
      // console.log('touching pos id',id);
  });
  var pos = 40;
  var length = 0;
  var between_residues = 18;
  var pos_bend = 0;
  var bend = 0;
  var distance_between_rows = 30;
  var bend_direction = -1*orientation;
  if (term=='N') {
    residues.reverse();
  }
  // console.log('residues',residues);
  bezier_cache = [];
  animation_delay = (1000/residues.length);
  // console.log("redraw half");
  $.each(residues, function(key,val) {
    if (key==0) {
      // grab the first anchor point
      if (term=='N') {
        x1 = parseInt($('#snake>#'+(val+1)).attr('cx'));
        y1 = parseInt($('#snake>#'+(val+1)).attr('cy'));
      } else {
        x1 = parseInt($('#snake>#'+(val-1)).attr('cx'));
        y1 = parseInt($('#snake>#'+(val-1)).attr('cy'));
      }
      x2 = x1-90*orientation;
      if (term=='N') {
        y2 = y_min
        bezierY = (y_min+y1)/2+60*orientation;
      } else {
        y2 = y_max;
        bezierY = (y_max+y1)/2+60*orientation;
      }
      bezierX = x1+60*orientation;


      // console.log('move residue ',val,'attach to id',key,pos,length,y_max,"Helix anchor",$('#'+(val+1)),(val+1),x1,y1);

      if (term=='C' && y2<bezierY) y2=bezierY;

      points = "M "+(x1)+" "+(y1)+" Q"+(bezierX)+" "+(bezierY)+" "+(x2)+" "+(y2);
      length = lengthbezier([x1,y1],[bezierX,bezierY],[x2,y2],0.001);

      // console.log("new path",points,"length",length);
      // $("path."+term+"-term.long").remove();
      // $("line."+term+"-term.long").remove();
      // path = "<path class='"+term+"-term long' d='" + points + "' stroke='black' fill='none' stroke-width='4' />";
      // $("#snake").append(path);
      $("rect."+term+"-term.long").attr('y',(y1+y2)/2);
      $("rect."+term+"-term.long").next().attr('y',13+(y1+y2)/2);

    }
    if (pos<length) {
      // on first bend
      [where,cur_length] = wherebezier([x1,y1],[bezierX,bezierY],[x2,y2],0.001,pos);
      // console.log(id,key,where,cur_length);

      if (key==0){
        $("line."+term+"-term.long").attr('x1',x1);
        $("line."+term+"-term.long").attr('y1',y1);
        $("line."+term+"-term.long").attr('x2',where[0]);
        $("line."+term+"-term.long").attr('y2',where[1]);
        // line = "<line class='"+term+"-term long' x1="+x1+" y1="+y1+" x2="+where[0]+" y2="+where[1]+" stroke='black' fill='none' stroke-width='2' />"
        // $("#snake").prepend(line);
      }

    } else {
        if (where[1]<y2 && term=='C') where[1]=y2;
        if (where[1]>y2 && term=='N') where[1]=y2;
        if (pos_bend==0 && bend!=0){
          where[0] = where[0]-between_residues*bend_direction;
          where[1] = where[1]+orientation*distance_between_rows/2;
        } else if (pos_bend==between_residues && bend!=0) {
          //#if 2nd residue in line put in middle
           where[0] = where[0]+between_residues*bend_direction;
           where[1] = where[1]+orientation*distance_between_rows/2;
         } else {
          where[0] = where[0]+between_residues*bend_direction;
          where[1] =  where[1];
        }
        last_bend_x = where[0];
        last_bend_y = where[1];
        pos_bend += between_residues;


        if (pos_bend>=Math.abs(x2-x_max)-40) {
          //no more bend left
            pos_bend = 0;
            bend += 1;
            if (bend_direction==1){
                bend_direction = -1;
            } else if (bend_direction==-1){
                bend_direction = 1;
            }
          }
    }
      // move residue
      if (where[0]!=$("#"+val).attr('cx') && where[1]!=$("#"+val).attr('cy')) {
        $("#"+val).attr('cx',where[0]);
        $("#"+val).attr('cy',where[1]);
        $("#"+val+"t").attr('x',where[0]);
        $("#"+val+"t").attr('y',where[1]+6);

        // delay = Math.round(key*animation_delay/10)*10;

        // $("#"+val).hide().delay( delay ).show(1);
        // $("#"+val+"t").hide().delay( delay ).show(1);

        if (where[1]<y_min && term=="N") y_min = where[1];
        if (where[1]>y_max && term=="C") y_max = where[1];
        // console.log('delay!',key*100,animation_delay,key*animation_delay,delay);
      }
      // $("#"+val).delay( key*100 ).queue(function() {$(this).attr('cx',where[0]).dequeue(); });
      // $("#"+val).delay( key*100 ).queue(function() { $(this).attr('y',where[1]).dequeue(); });
      // $("#"+val+"t").delay( key*100 ).queue(function() { $(this).attr('x',where[0]).dequeue(); });
      // $("#"+val+"t").delay( key*100 ).queue(function() { $(this).attr('y',where[1]+6).dequeue(); });

    pos += between_residues;

  });
  // console.log("redraw done");
      // move label box
  // console.log(term,"min",y_min,"max",y_max);

}

function applyPresentColors(target) {

    // Color all residues by their amino acid using the presetColors array

    $('#'+target).find(".rcircle").each(function( index ){
          aa =  $(this).next().text().trim();
          $(this).css("fill", presetColors[aa][0]);
          $(this).next().css("fill", presetColors[aa][1]);
        });

};

function resetColors(target) {

    // Reset color of all residues

    $('#'+target).find(".rcircle").each(function( index ){
          aa =  $(this).next().text();
          $(this).css("fill", 'white');
          $(this).next().css("fill", 'black');
        });

}

function maxmin() {
    // console.log('maxmin start');
    margin = 50;
    svgmax = 0;
    svgmin = 0;
    x_svgmax = 0;
    count = 0;
    classmax = '';
    classmin = '';
    counter = 0;
    if (!$('#snake').length) return
    // console.log("temp",y_max,y_min);
    $('#snake').children('.rtext').each(function () {
        counter += 1;
        y = parseInt($(this).attr( "y" ));
        x = parseInt($(this).attr( "x" ));
        classtext = $(this).attr( "class" );
        // test = $(this).attr("original_title");
        // test2 = $(this).css("display");
        if ($(this).css("display")!='none') {
            count = count +1;
            if (y<svgmin) {
                svgmin = y;
                classmin = classtext;
                }
            if (y>svgmax) {
                classmax = classtext;
                svgmax= y;
             }
            if (x>x_svgmax) x_svgmax = x;

        }
    });

    // Alternative for height - otherwise the height of the element itself is skipped
    //svgmax = document.getElementById('snake').getBBox().height;

    // if (svgmin>y_min) svgmin = y_min;
    // if (svgmax<y_max) svgmax = y_max;

    // console.log('max '+svgmax+' '+classmax+' min'+svgmin+' '+classmin+' count'+count);

    var svg = $('#snake').closest('svg')[0];
    check = document.getElementById("snakeplot").getAttribute('viewBox');
    // console.log('check',check)
    if (typeof check !== typeof undefined && check !== false && check !== null ) {
      oldheight = check.split(" ")[3];
      width = check.split(" ")[2];
    } else {
      // console.log('not found it');
      oldheight = $(svg).attr('height');
      width = $("#snakeplot").attr('width');
    }

    newheight = (svgmax-svgmin+margin*2)+margin/2; // add extra border/padding

    // console.log('New height:'+ newheight +' old height:'+oldheight);
    // console.log("Prev attr"+$('#snake').attr("transform"));
    if (newheight!=oldheight) {
        svg.setAttribute('height', (svgmax-svgmin+margin*2));
        $('#snake').attr("transform", "translate(0," + (-svgmin+margin) + ")");

        // $('#snakeplot')[0].attr("viewBox", "0 0 " + width + " " + newheight);
        document.getElementById("snakeplot").setAttribute("viewBox", "0 0 " + width + " " + newheight);


        svg.setAttribute('height', "100%");
        svg.setAttribute('width', "100%");
    }

    // console.log("New attr"+$('#snake').attr("transform"));
    // console.log('maxmin done');
}

$( document ).ready(function() {
    // var elements = document.getElementsByClassName('long')

    // for (var i = 0; i < elements.length; i++){
    //     elements[i].style.display = 'none';
    // }
    y_min = 0;
    y_max = 0;
    // maxmin();
    $(".long").hide();


    $('rect').each(function(){

        rectclass = $(this).attr('class');
        if (rectclass) {
            if (rectclass.indexOf("CL") >= 0 && rectclass.indexOf("long") >= 0) {

                numResidues = ($('.'+rectclass.replace(/ /g,".")).length-3)/2

                // console.log('class:'+rectclass+' count'+numResidues);

                if (numResidues<10) {
                    toggleLoop('.'+rectclass.split(' ')[0],'',1);
                }
            }
        }
    });

    maxmin();

    $("text").tooltip({
        'container': 'body',
        'placement': 'top',
        'animation': false,
        'html' : true
    });


    $("circle").tooltip({
        'container': 'body',
        'placement': 'top',
        'animation': false,
        'html' : true
    });

    $("circle").hover(function(){
        $('.tooltip').css('top',parseInt($('.tooltip').css('top')) + 2.8 + 'px')
    });

});

$(".rtext").click(function() {
    parentid = $(this).closest('svg').attr('id');
    newcolor = $(".pick-color."+parentid+".selected").attr('id');
    if (newcolor) {
      newcolor = newcolor.split('-');
    } else {
      custom =  $("#custom_color_"+parentid).val();
      custom_text = getContrast50(custom);
      newcolor = ["pick",custom,custom_text];
    }
  console.log(newcolor);

  $(this).css("fill", newcolor[2]);
  $(this).prev().css("fill", newcolor[1]);
});

$(".rcircle").click(function() {
    parentid = $(this).closest('svg').attr('id');
    newcolor = $(".pick-color."+parentid+".selected").attr('id');

    if (newcolor) {
      newcolor = newcolor.split('-');
    } else {
      custom =  $("#custom_color_"+parentid).val();
      custom_text = getContrast50(custom);
      newcolor = ["pick",custom,custom_text];
    }
  console.log(newcolor);
  $(this).css("fill", newcolor[1]);
  $(this).next().css("fill", newcolor[2]);
});



$("#snake_svg_link").click(function() {
    svgAsDataUri(document.getElementById("snakeplot"),{}, function(uri) {
        $("#snake_svg_link").attr('href',uri);
    });
});
$("#snake_svg_link2").click(function() {
    svgAsDataUri(document.getElementById("snakeplot2"),{}, function(uri) {
        $("#snake_svg_link2").attr('href',uri);
    });
});
$("#helix_svg_link").click(function() {
    svgAsDataUri(document.getElementById("helixbox"),{}, function(uri) {
        $("#helix_svg_link").attr('href',uri);
    });
});

function ajaxMutants(plotid,protein) {

  resetColors(plotid);

    $.getJSON( '/mutations/ajax/'+protein+'/', function( data ) {
      $.each( data, function( key, val ) {

         var ligands = [], bigincreases=0, increases = 0, bigdecreases=0, decreases = 0, unchanged=0, unspecified = 0;


         $.each( val, function( key, v ) {
          if( !(ligands[v[1]]) ) ligands[v[1]] = [];
          ligands[v[1]].push(v[0])
          if (v[0]>10) {
            bigincreases ++; //mix-up increase is decrease.
          } else if (v[0]>5) {
            increases ++;
          } else if (v[0]>0) {
            unchanged ++;
          }  else if (v[0]<-10) {
            bigdecreases ++;
          } else if (v[0]<-5) {
            decreases ++;
          } else if (v[0]<0) {
            unchanged ++;
          } else if (v[2]=='No effect on') {
            unchanged ++;
          } else if (v[2]=='No effect') {
            unchanged ++;
          } else if (v[2]=='Abolish') {
            bigincreases ++;
          } else if (v[2]=='Abolished effect') {
            bigincreases ++;
          } else if (v[2]=='Gain of') {
            bigdecreases ++;
          } else if (v[2]=='Increase') {
            increases ++;
          } else if (v[2]=='Decrease') {
            decreases ++;
          } else {
            unspecified ++;
          }
         });

         extra = "\n" + String(val[0].length) + " mutations: " +
          (decreases+bigdecreases) +" increases | " +
         (increases+bigincreases) +" decreases  |  " +
          (unchanged) +" Unchanged | " +
          unspecified + " Unspecified";

        counts = [(increases+bigincreases),(decreases+bigdecreases),(unchanged)];
        winner = counts.indexOf(Math.max.apply(window,counts));
        winner2 = Math.max.apply(window,counts);
        color = "#D9D7CE";
        color_letter = "#000";
        // if (winner==0 && winner2) {
        //   if (increases>bigincreases) {
        //     color = "#FF7373";
        //     color_letter = "#FFF";
        //   } else {
        //     color = "#FA1111";
        //     color_letter = "#FDFF7B";
        //   }
        // } else if (winner==1) {
        //   if (decreases>bigdecreases) {
        //     color = "#87E88F";
        //   } else {
        //     color = "#66B36C";
        //   }
        // } else if (winner==2) {
        //   color = "#F7DA00";
        //   color_letter = "#000";
        // }

        if (bigincreases>0) {
            color = "#FF7373";
            color_letter = "#FFF";
        } else if (increases>0) {
            color = "#FA1111";
            color_letter = "#FDFF7B";
        } else if (bigdecreases>0) {
            color = "#66B36C";
        } else if (decreases>0) {
            color = "#87E88F";
        } else  {
          color = "#F7DA00";
          color_letter = "#000";
        }

        console.log(counts + " " + counts.indexOf(Math.max.apply(window,counts)));


         original_title = $('#'+plotid).find("#"+key).attr('original_title')

         $('#'+plotid).find("#"+key).css("fill", color);
         $('#'+plotid).find("#"+key).next().css("fill",color_letter );
         $('#'+plotid).find("#"+key).attr('title',original_title+extra);
         $('#'+plotid).find("#"+key+"t").attr('title',original_title+extra);


      });
    $("circle").tooltip('fixTitle');
    $("text").tooltip('fixTitle');

    });
}

function ajaxMutantsPos(plotid) {

  resetColors(plotid);

    var pos = jQuery.parseJSON(mutant_json);

    $.each( pos, function( key, val ) {
         var ligands = [], bigincreases=0, increases = 0, bigdecreases=0, decreases = 0, unchanged=0, unspecified = 0;


         $.each( val[0], function( key, v ) {
          if( !(ligands[v[1]]) ) ligands[v[1]] = [];
          ligands[v[1]].push(v[0])
          if (v[0]>10) {
            bigincreases ++; //mix-up increase is decrease.
          } else if (v[0]>5) {
            increases ++;
          } else if (v[0]>0) {
            unchanged ++;
          } else if (v[0]<-10) {
            bigdecreases ++;
          } else if (v[0]<-5) {
            decreases ++;
          } else if (v[0]<0) {
            unchanged ++;
          } else if (v[2]=='No effect on') {
            unchanged ++;
          } else if (v[2]=='No effect') {
            unchanged ++;
          } else if (v[2]=='Abolish') {
            bigincreases ++;
          } else if (v[2]=='Abolished effect') {
            bigincreases ++;
          } else if (v[2]=='Gain of') {
            bigdecreases ++;
          } else if (v[2]=='Increase') {
            increases ++;
          } else if (v[2]=='Decrease') {
            decreases ++;
          } else {
            unspecified ++;
          }
         });

         extra = "\n" + String(val[0].length) + " mutations: " +
          (decreases+bigdecreases) +" increases | " +
         (increases+bigincreases) +" decreases  |  " +
          (unchanged) +" Unchanged | " +
          unspecified + " Unspecified";

        counts = [(increases+bigincreases),(decreases+bigdecreases),(unchanged)];
        winner = counts.indexOf(Math.max.apply(window,counts));
        winner2 = Math.max.apply(window,counts);
        color = "#D9D7CE";
        color_letter = "#000";

        if (bigincreases>0) {
          color = "#FA1111";
          color_letter = "#FDFF7B";
        } else if (increases>0) {
          color = "#FF7373";
          color_letter = "#FFF";
        } else if (bigdecreases>0) {
            color = "#66B36C";
        } else if (decreases>0) {
            color = "#87E88F";
        } else  {
          color = "#F7DA00";
          color_letter = "#000";
        }

         original_title = $('#'+plotid).find("#"+key).attr('original_title')
         $('#'+plotid).find("#"+key).css("fill", color);
         $('#'+plotid).find("#"+key).next().css("fill",color_letter );
         $('#'+plotid).find("#"+key).attr('title',original_title+extra);
         $('#'+plotid).find("#"+key+"t").attr('title',original_title+extra);


      });
    $("circle").tooltip('fixTitle');
    $("text").tooltip('fixTitle');


}

function ajaxInteractions(plotid,protein) {

  resetColors(plotid);

    $.getJSON( '/interaction/ajax/'+protein+'/', function( data ) {
      $.each( data, function( key, val ) {

        var flags = [], falgsAA = [], output = [], outputAA = [], l = val.length, i;
        for( i=0; i<l; i++) {
            if( flags[val[i][1]]) continue;
            flags[val[i][1]] = true;
            output.push(val[i][1]);
        }
        for( i=0; i<l; i++) {
            if( flags[val[i][0]]) continue;
            flags[val[i][0]] = true;
            outputAA.push(val[i][0]);
        }

         extra = "\n" + String(val.length) + " interactions | Type: "+ output +" | Residue in crystal:"+ outputAA;


         $('#'+plotid).find("#"+key).css("fill", "#E60A0A");
         $('#'+plotid).find("#"+key).next().css("fill", "#FDFF7B");

         original_title = $('#'+plotid).find("#"+key).attr('original_title')


         $('#'+plotid).find("#"+key).attr('title',original_title+extra);
         $('#'+plotid).find("#"+key+"t").attr('title',original_title+extra);


      });
    $("circle").tooltip('fixTitle');
    $("text").tooltip('fixTitle');

    });
}

function ajaxBarcode(plotid,protein) {

  resetColors(plotid);

    var textboxvalue = $('input[name=cutoff]').val();
    $.getJSON( '/signprot/ajax/barcode/'+protein+'/'+textboxvalue, function( data ) {
      $.each( data, function( key, val ) {

        if (val[1]=='Conserved') {
            color = "#b162a7";
            color_letter = "#fefdfd";
            $('#'+plotid).find("#"+key).next().css("fill", color_letter);
            extra = "\n" + String(val[1]);
        } else if (val[1]=='Evolutionary neutral') {
            color = "#f8dfb4";
            extra = "\n" + String(val[1]);
        } else if (val[1]=='NA') {
            color = "#ffffff"
            extra = "\n" + String(val[1]);
        } else  {
            color = "#4dc7e6";
            extra = "\n" + String(val[1]);
        }

         original_title = $('#'+plotid).find("#"+key).attr('original_title')

         $('#'+plotid).find("#"+key).css("fill", color);
         $('#'+plotid).find("#"+key).attr('title',original_title+extra);
         $('#'+plotid).find("#"+key+"t").attr('title',original_title+extra);

      });
    $("circle").tooltip('fixTitle');
    $("text").tooltip('fixTitle');

    });
}
function ajaxCancerMutation(plotid, protein) {

  resetColors(plotid);

    $.getJSON( '/mutational_landscape/ajax/CancerMutation/'+protein+'/', function( data ) {
      $.each( data, function( key, val ) {
        // NM.allele_frequency, NM.allele_count, NM.allele_number, NM.number_homozygotes
         extra = "\nAAchange: " + "-->" + String(val[0]);

         color = "#7572b1";
         color_letter = "#fefdfd";
         $('#'+plotid).find("#"+key).next().css("fill", color_letter);

         original_title = $('#'+plotid).find("#"+key).attr('original_title')
         $('#'+plotid).find("#"+key).css("fill", color);
         $('#'+plotid).find("#"+key).attr('title',original_title+extra);
         $('#'+plotid).find("#"+key+"t").attr('title',original_title+extra);


      });
    $("circle").tooltip('fixTitle');
    $("text").tooltip('fixTitle');

    });
}

function ajaxDiseaseMutation(plotid, protein) {

  resetColors(plotid);

    $.getJSON( '/mutational_landscape/ajax/DiseaseMutation/'+protein+'/', function( data ) {
      $.each( data, function( key, val ) {
        // NM.allele_frequency, NM.allele_count, NM.allele_number, NM.number_homozygotes
         extra = "\nAAchange: " + "-->" + String(val[0]);

         color = "#52133b";
         color_letter = "#fefdfd";
         $('#'+plotid).find("#"+key).next().css("fill", color_letter);

         original_title = $('#'+plotid).find("#"+key).attr('original_title')
         $('#'+plotid).find("#"+key).css("fill", color);
         $('#'+plotid).find("#"+key).attr('title',original_title+extra);
         $('#'+plotid).find("#"+key+"t").attr('title',original_title+extra);


      });
    $("circle").tooltip('fixTitle');
    $("text").tooltip('fixTitle');

    });
}
function ajaxNaturalMutation(plotid, protein) {

  resetColors(plotid);

    $.getJSON( '/mutational_landscape/ajax/NaturalMutation/'+protein+'/', function( data ) {
      $.each( data, function( key, val ) {
         extra = "\nType: " + "-->" + String(val[5]) +
         "\nAAchange: " + "-->" + String(val[0]) +
         "\nAllele Frequency: " + String(val[1]) +
         "\nAllele Count: " + String(val[2]) +
         "\nAllele Number: " + String(val[3]) +
        "\nNumber of Homozygotes: " + String(val[4]) +
        "\nFunctional Annotation: " + String(val[8]) +
        "\nPredicted effect (SIFT/PolyPhen): <span style='color:"+String(val[7])+"'> "+ String(val[6]);

        //  color = "#c40100";
         $('#'+plotid).find("#"+key).next().css("fill", "#fefdfd");
         original_title = $('#'+plotid).find("#"+key).attr('original_title')
         $('#'+plotid).find("#"+key).css("fill", String(val[7]));
         $('#'+plotid).find("#"+key).attr('title',original_title+extra);
         $('#'+plotid).find("#"+key+"t").attr('title',original_title+extra);


      });
    $("circle").tooltip('fixTitle');
    $("text").tooltip('fixTitle');

    });
}

function ajaxNaturalMutationPos(plotid) {

  resetColors(plotid);

    var pos = jQuery.parseJSON(natural_mutations_json);

    var color_code = pos['color']

      $.each(pos, function( key, val ) {
         //console.log("Yes", pos);

         extra = "\nVariants: " + "-->" + String(val['AA']) +
        "\nNumber of Proteins: " + String(val['val']);

         if (val['val']==0) {
            color_letter = "#000000"
         } else  {
            color_letter = "#fefdfd";
         }

         color = color_code[val['val']-1];
         $('#'+plotid).find("#"+key).next().css("fill", color_letter);

         original_title = $('#'+plotid).find("#"+key).attr('original_title')
         $('#'+plotid).find("#"+key).css("fill", color);
         $('#'+plotid).find("#"+key).attr('title',original_title+extra);
         $('#'+plotid).find("#"+key+"t").attr('title',original_title+extra);

      });
    $("circle").tooltip('fixTitle');
    $("text").tooltip('fixTitle');
}

function ajaxPTMs(plotid, protein) {

  resetColors(plotid);

    $.getJSON( '/mutational_landscape/ajax/PTM/'+protein+'/', function( data ) {
      $.each( data, function( key, val ) {
         extra = "\nModification: " + String(val[0]);

         $('#'+plotid).find("#"+key).next().css("fill", "#fefdfd");
         original_title = $('#'+plotid).find("#"+key).attr('original_title')
         $('#'+plotid).find("#"+key).css("fill", "#000000");
         $('#'+plotid).find("#"+key).attr('title',original_title+extra);
         $('#'+plotid).find("#"+key+"t").attr('title',original_title+extra);


      });
    $("circle").tooltip('fixTitle');
    $("text").tooltip('fixTitle');

    });
}

function ajaxPTMPos(plotid) {

  resetColors(plotid);

    var pos = jQuery.parseJSON(ptms_json);

    var color_code = pos['color']

      $.each(pos, function( key, val ) {

         extra = "\nModifications: " + String(val['mod']);

         $('#'+plotid).find("#"+key).next().css("fill", "#fefdfd");

         original_title = $('#'+plotid).find("#"+key).attr('original_title')
         $('#'+plotid).find("#"+key).css("fill", "#000000");
         $('#'+plotid).find("#"+key).attr('title',original_title+extra);
         $('#'+plotid).find("#"+key+"t").attr('title',original_title+extra);

      });
    $("circle").tooltip('fixTitle');
    $("text").tooltip('fixTitle');
}

function ajaxCancerMutationPos(plotid) {

  resetColors(plotid);

    var pos = jQuery.parseJSON(cancer_mutations_json);

    var color_code = pos['color']

      $.each(pos, function( key, val ) {

         extra = "\nAAchanges: " + "-->" + String(val['AA']) +
        "\nNumber of Proteins: " + String(val['val']);

         if (val['val']==0) {
            color_letter = "#000000"
         } else  {
            color_letter = "#fefdfd";
         }

         color = color_code[val['val']-1];
         $('#'+plotid).find("#"+key).next().css("fill", color_letter);

         original_title = $('#'+plotid).find("#"+key).attr('original_title')
         $('#'+plotid).find("#"+key).css("fill", color);
         $('#'+plotid).find("#"+key).attr('title',original_title+extra);
         $('#'+plotid).find("#"+key+"t").attr('title',original_title+extra);


      });
    $("circle").tooltip('fixTitle');
    $("text").tooltip('fixTitle');

}

function ajaxDiseaseMutationPos(plotid) {

  resetColors(plotid);

    var pos = jQuery.parseJSON(disease_mutations_json);

    var color_code = pos['color']

      $.each(pos, function( key, val ) {

         extra = "\nAAchanges: " + "-->" + String(val['AA']) +
        "\nNumber of Proteins: " + String(val['val']);

         color_letter = "#fefdfd";
         color = color_code[val['val']-1];
         $('#'+plotid).find("#"+key).next().css("fill", color_letter);

         original_title = $('#'+plotid).find("#"+key).attr('original_title')
         $('#'+plotid).find("#"+key).css("fill", color);
         $('#'+plotid).find("#"+key).attr('title',original_title+extra);
         $('#'+plotid).find("#"+key+"t").attr('title',original_title+extra);


      });
    $("circle").tooltip('fixTitle');
    $("text").tooltip('fixTitle');
}

function ajaxInterface(plotid,protein) {

  resetColors(plotid);

    $.getJSON( '/signprot/ajax/interface/'+protein+'/', function( data ) {
      $.each( data, function( key, val ) {

         extra = "\n" + String("Receptor interface position");
         $('#'+plotid).find("#"+key).css("fill", "#e60a0a");
         $('#'+plotid).find("#"+key).next().css("fill", "#fafb74");

         original_title = $('#'+plotid).find("#"+key).attr('original_title')

         $('#'+plotid).find("#"+key).attr('title',original_title+extra);
         $('#'+plotid).find("#"+key+"t").attr('title',original_title+extra);


      });
    $("circle").tooltip('fixTitle');
    $("text").tooltip('fixTitle');

    });
}

function ajaxInteractionsPos(plotid) {

  resetColors(plotid);


  var pos = jQuery.parseJSON(interaction_json);

      $.each( pos, function( key, val ) {

        var flags = [], falgsAA = [], output = [], outputAA = [], l = val.length, i;
        for( i=0; i<l; i++) {
            if( flags[val[i][1]]) continue;
            flags[val[i][1]] = true;
            output.push(val[i][1]);
        }
        for( i=0; i<l; i++) {
            if( flags[val[i][0]]) continue;
            flags[val[i][0]] = true;
            outputAA.push(val[i][0]);
        }

         extra = "\n" + String(val.length) + " interactions | Type: "+ output +" | Residue in crystal:"+ outputAA;


         $('#'+plotid).find("#"+key).css("fill", "#E60A0A");
         $('#'+plotid).find("#"+key).next().css("fill", "#FDFF7B");

         original_title = $('#'+plotid).find("#"+key).attr('original_title')


         $('#'+plotid).find("#"+key).attr('title',original_title+extra);
         $('#'+plotid).find("#"+key+"t").attr('title',original_title+extra);


      });
    $("circle").tooltip('fixTitle');
    $("text").tooltip('fixTitle');

}

function construct_annotations(plotid) {

  resetColors(plotid);

  var pos = jQuery.parseJSON(annotations_json);

      $.each( pos, function( key, val ) {

        // var flags = [], falgsAA = [], output = [], outputAA = [], l = val.length, i;
        // for( i=0; i<l; i++) {
        //     if( flags[val[i][1]]) continue;
        //     flags[val[i][1]] = true;
        //     output.push(val[i][1]);
        // }
        // for( i=0; i<l; i++) {
        //     if( flags[val[i][0]]) continue;
        //     flags[val[i][0]] = true;
        //     outputAA.push(val[i][0]);
        // }

        //  extra = "\n" + String(val.length) + " interactions | Type: "+ output +" | Residue in crystal:"+ outputAA;
         extra = "<br>"+val[1]; //.replace(/<br>/g, '&#013;');


         $('#'+plotid).find("#"+key).css("fill", val[2]);
         $('#'+plotid).find("#"+key).next().css("fill", val[3]);

         original_title = $('#'+plotid).find("#"+key).attr('original_title')

         $('#'+plotid).find("#"+key).attr('title',original_title+extra);
         $('#'+plotid).find("#"+key+"t").attr('title',original_title+extra);


      });
    $("circle").tooltip('fixTitle');
    $("text").tooltip('fixTitle');

}

function ajaxInteractionsLigand(protein,ligand) {

    resetColors('snakeplot');
    resetColors('helixbox');

    $.getJSON( '/interaction/ajaxLigand/'+protein+'/'+ligand, function( data ) {
      $.each( data, function( key, val ) {

        var flags = [], falgsAA = [], output = [], outputAA = [], l = val.length, i;
        for( i=0; i<l; i++) {
            if( flags[val[i][1]]) continue;
            flags[val[i][1]] = true;
            output.push(val[i][1]);
        }
        for( i=0; i<l; i++) {
            if( flags[val[i][0]]) continue;
            flags[val[i][0]] = true;
            outputAA.push(val[i][0]);
        }

         extra = "\n" + String(val.length) + " interactions | Type: "+ output +" | Residue in crystal:"+ outputAA;


         $('[id='+key+']').css("fill", "#E60A0A");
         $('[id='+key+']').next().css("fill", "#FDFF7B");

         original_title = $("#"+key).attr('original_title')

         $('[id='+key+']').attr('title',original_title+extra);
         $('[id='+key+'t]').attr('title',original_title+extra);


      });
    $("circle").tooltip('fixTitle');
    $("text").tooltip('fixTitle');

    });
}


$(".pick-color").click(function() {
    plottype = $(this).attr('class').split(' ')[1];

    console.log($(this).attr('id'));
    $(".pick-color."+plottype).css('borderWidth','2px');
    $(".pick-color."+plottype).css('height','20px');
    $(".pick-color."+plottype).removeClass('selected');
    $(this).css('borderWidth','3px');
    $(this).css('height','22px');
    $(this).addClass('selected');

});

$( document ).ready(function() {
  if ( $( "#cp2_helixbox" ).length ) {
    $('#cp2_helixbox').colorpicker();
  }
  if ( $( "#cp2_snakeplot" ).length ) {
    $('#cp2_snakeplot').colorpicker();
  }
});

function getContrast50(hexcolor){
     return (parseInt(hexcolor.replace('#', ''), 16) > 0xffffff/2) ? 'black':'white';
}

function reload_tooltips() {
    $("text").tooltip({
        'container': 'body',
        'placement': 'top',
        'animation': false,
        'html' : true
    });


    $("circle").tooltip({
        'container': 'body',
        'placement': 'top',
        'animation': false,
        'html' : true
    });

    $("circle").hover(function(){
        $('.tooltip').css('top',parseInt($('.tooltip').css('top')) + 2.8 + 'px')
    });
}
