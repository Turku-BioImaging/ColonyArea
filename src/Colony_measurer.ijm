/** Copyright 2013, Camilo Guzmán, Manish Bagga and Daniel Abankwa
    

    This file is part of ColonyArea.

    ColonyArea is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ColonyArea is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ColonyArea.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
This source code is maintained by Turku BioImaging
https://github.com/Turku-BioImaging/ColonyArea
https://bioimaging.fi
*/


print("This window will be cleared!!");
selectWindow("Log");
run("Close");
print("Well #         Area Percent        Intensity Percent");
run("Set Measurements...", "area mean area_fraction redirect=None decimal=1");
dir=getInfo("image.directory");
name=getInfo("image.filename");
height=getHeight();
width=getWidth();
in=0;
out=0;
getSelectionBounds(x,y,h,w);

//print("height"+height+"   width"+width+"   x"+x+"   y"+y+"   h"+h+"   w"+w);
in=0;
out=0;
total=0;
total1=0;
dist = 0
for (i=x;i<x+h;i++)
	{for (j=y;j<y+w;j++)
		{total1=total1+1;
		dist = (i-height/2)*(i-height/2)+(j-width/2)*(j-width/2);
//print(dist);	
	if ((i-height/2)*(i-height/2)*width*width/4+(j-width/2)*(j-width/2)*height*height/4<=(height*height*width*width/16))
		{in=in+1;}
		else
		{out=out+1;}


		}
	}
total=h*w;
ratio=in/total;
//print("  in"+in+" out"+out+" ratio"+ratio);
from=1;
to=nSlices;
r=0;
if(h*w<height*width) 
{r=1;
from=getNumber("Enter the well # from which to measure: ", getSliceNumber());
to=getNumber("Enter the well # till which to measure: ", getSliceNumber());
}
if(from<1||to<1||from>nSlices||to>nSlices||from>to)
{print("Please enter the limits correctly!");
return;}
for (n=from; n<=to; n++)

{setSlice(n);run("Measure");
areafrac=getResult("%Area");
area=getResult("Area");
mean=getResult("Mean");
areafrac=areafrac/ratio;

sumofi=area*mean;
maxsumofi=area*255*ratio;
intensityfrac=sumofi*100/maxsumofi;
if(areafrac>=9.995)
{print(n+"                  "+d2s(areafrac,2)+"                    "+d2s(intensityfrac,2));}
else
{print(n+"                  "+d2s(areafrac,2)+"                      "+d2s(intensityfrac,2));}
}
setSlice(1);
selectWindow("Log");
if(r==1)
{print(" \n"+"x="+x+" y="+y+" h="+h+" w="+w);
}
name=replace(name, "\\.tif", "");
fullname=dir+"results_"+name+".txt";
i=0;
while(File.exists(fullname)==1)
{i=i+1;
fullname=dir+"results_"+name+"_"+i+".txt";}

//fullname=dir+"results_"+name+".txt";
//print(fullname);
saveAs("Text", fullname);
selectWindow("Results");
run("Close");
