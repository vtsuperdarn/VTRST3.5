/* hlpstr.h
   ========
   Author: R.J.Barnes
*/

/*
 $License$
*/


char *hlpstr[]={
"fov_plot - Plot the radar fields of view.\n",
"fov_plot --help\n",
"fov_plot [options]\n",
"fov_plot -ps [-xp xoff] [-yp yoff] [options]\n",
"fov_plot -ppm [options]\n",
"fov_plot -ppmx [options]\n",
"fov_plot -xml [options]\n",
"fov_plot -x [-display display] [-xoff xoff] [-yoff yoff] [options]\n",

"--help\tprint the help message and exit.\n",
"-cf cfname\tread command line options from the file cfname.\n",
"-wdt width\tset the width of the plot to width.\n",
"-hgt height\tset the height of the plot to height.\n",
"-pad pad\tset the paddding around the edges of the plot to pad.\n",
"-sf scale\tset the scale factor (magnification) to scale. The default scale factor is 1.\n",
"-ortho\tuse an orthographic projection.\n",
"-stereo\tuse a stereographic projection.\n",
"-gvp\tuse a general vertical projection.\n",
"-lat lat\tcenter the plot on the latitude given by lat.\n",
"-lon lon\tcenter the plot on the longitude given by lon.\n",
"-latmin latmin\tadjust the scale factor so that the lowest visible latitude is  latmin. Applies when the stereographic projection is used.\n",
"-mag\tuse magnetic coordinates.\n",
"-rotate\trotate the plot so that the local noon is at the top of the plot.\n",
"-flip\tflip the direction of the X-axis.\n",
"-square\tforce the use of a square bounding box around the plot.\n",
"-coast\tplot coastlines.\n",
"-fcoast\tplot filled coastlines.\n",
"-bnd\tplot state boundaries.\n",
"-grd\tplot a grid.\n",
"-tmk\tplot a clock-dial grid showing the time.\n",
"-grdontop\tplot the grid on top of the image.\n",
"-fov\tplot the radar field of view.\n",
"-ffov\tplot the filled radar field of view.\n",
"-ofov\tplot the field of view of all the other radars.\n",
"-fofov\tplot the filled field of view of all the other radars.\n",
"-cfov\tplot the field of view of radars under construction.\n",
"-fcfov\tplot the filled field of view of radars under construction.\n",
"-dot\tMark the location of the radar site with a dot.\n",
"-dotr radius\tset the radius of the dot used to mark the station.\n",
"-tmtick tick\tset the grid interval for the time clock-dial to tick hours.\n",
"-lst\tuse local solar time rather than local time.\n",
"-term\tplot the terminator.\n",
"-term\tplot a shaded terminator.\n",
"-tmlbl\tLabel the time clock-dial.\n",
"-fontname fontname\tplot any labels using the font fontname.\n",
"-fontsize fontsize\tset the font size to fontsize points.\n",
"-lnewdt lnewdt\tset the line width to lnewdt.\n",
"-bgcol rrggbb\tset the background color to rrggbb, specified as the hexadecimal value for the 24-bit red,green and blue component color.\n",
"-txtcol rrggbb\tset the color of the text to rrggbb, specified as the hexadecimal value for the 24-bit red,green and blue component color.\n",
"-grdcol rrggbb\tset the color of the grid to rrggbb, specified as the hexadecimal value for the 24-bit red,green and blue component color.\n",
"-cstcol rrggbb\tset the color of the coastline to rrggbb, specified as the hexadecimal value for the 24-bit red,green and blue component color.\n",
"-bndcol rrggbb\tset the color of the state boundaries to rrggbb, specified as the hexadecimal value for the 24-bit red,green and blue component color.\n",
"-lndcol rrggbb\tset the color of the land to rrggbb, specified as the  hexadecimal value for the 24-bit red,green and blue component color.\n",
"-seacol rrggbb\tset the color of the sea to rrggbb, specified as the  hexadecimal value for the 24-bit red,green and blue component color.\n",
"-trmcol rrggbb\tset the color of the terminator outline to rrggbb, specified as the hexadecimal value for the 24-bit red,green and blue component color.\n",
"-ftrmcol rrggbb\tset the color of the filled terminator to rrggbb, specified as the hexadecimal value for the 24-bit red,green and blue component color.\n",
"-tmkcol rrggbb\tset the color of the time clock-dial to rrggbb, specified as the hexadecimal value for the 24-bit red,green and blue component color.\n",
"-fovcol rrggbb\tset the color of the field of view outline to rrggbb, specified as the hexadecimal value for the 24-bit red,green and blue component color.\n",
"-ffovcol rrggbb\tset the color of the filled field of view to rrggbb, specified as the hexadecimal value for the 24-bit red,green and blue component color.\n",
"-ofovcol rrggbb\tset the color of the field of view outline of other radars to rrggbb, specified as the hexadecimal value for the 24-bit red,green and blue component color.\n",
"-fofovcol rrggbb\tset the color of the filled field of view of other radars to rrggbb, specified as the hexadecimal value for the 24-bit red,green and blue component color.\n",
"-cfovcol rrggbb\tset the color of the field of view outline of radars under construction to rrggbb, specified as the hexadecimal value for the 24-bit red,green and blue component color.\n",
"-fcfovcol rrggbb\tset the color of the filled field of view of radars under construction to rrggbb, specified as the hexadecimal value for the 24-bit red,green and blue component color.\n",
"-d yyyymmdd\tplot for the date yyyymmdd.\n",
"-t hr:mn\tplot for the time hr:mn.\n",
"-st stid\tplot the field of view for the radar with the station identifier stid.\n",
"-ps\tproduce a PostScript plot as the output.\n",
"-xp xoff\tset the X offset of the PostScript plot to xoff.\n",
"-yp xoff\tset the Y offset of the PostScript plot to yoff.\n",
"-ppm\tproduce a Portable PixMap (PPM) image as the output.\n",
"-ppmx\tproduce an extended Portable PixMap (PPMX) image as the output.\n",
"-xml\tproduce an XML image representation as the output.\n",
"-x\tplot on an X-terminal.\n",
"-display display\tconnect to the xterminal named display.\n",
"-xoff xoff\topen the window, xoff pixels from the left edge of the screen.\n",
"-yoff yoff\topen the window ypad pixels from the top edge of the screen.\n",

NULL};