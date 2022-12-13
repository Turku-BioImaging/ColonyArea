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


  macro "Manual colony thresholder" {
print("This window will be cleared.");
selectWindow("Log");
run("Close");
       dir=getInfo("image.directory");
dir1=dir;       
name=getInfo("image.filename");
name1=name;

       run("8-bit");
       run("Set Measurements...", "area mean area_fraction  redirect=None decimal=1"); 
t1=newArray(39);
t2=newArray(39);
t3=newArray(39);
t4=newArray(39);
t5=newArray(39);
s1=newArray(37);
s2=newArray(37);
s3=newArray(37);
s4=newArray(37);
s5=newArray(37);
       final=newArray(nSlices);
	intensity=newArray(nSlices);
      again=0;
avg=0; 
denom=0;
iteration=1;
default_w=getNumber("enter reference well number:",getSliceNumber());
setSlice(default_w);
default=0;
for (upper=30;upper<=220;upper=upper+1)
{setThreshold(1,upper);
run("Measure");
if(getResult("%Area")<=78){default=upper;}
}
resetThreshold();

default_t=getNumber("enter reference threshold:",default-40);
max_i=newArray(nSlices);
threshold=newArray(nSlices);
	print("Well_#     Threshold");
//                                            run("Image Sequence... ", "format=TIFF name=["+name+"] start=1 digits=1 save=["+dir+"]");
       
for (n=1; n<=nSlices; n++) {
		selectWindow(name1);
          setSlice(n);
          fullname="Well "+n+" of "+name;
          m=0;
	max_intensity=0;
min_intensity=20;
//if(iteration==1) {         
for (upper=30; upper<=220; upper=upper+5) {
             setThreshold(1, upper);       
             run("Measure");
             t1[m]=getResult("%Area");
		setThreshold(1,upper+1);
		run("Measure");
		t2[m]=getResult("%Area");
		setThreshold(1,upper+2);
		run("Measure");
		t3[m]=getResult("%Area");
		setThreshold(1,upper+3);
		run("Measure");
		t4[m]=getResult("%Area");
		setThreshold(1,upper+4);
		run("Measure");
		t5[m]=getResult("%Area");
             
if(t1[m]<=0.7)
{
min_intensity=upper;

}


if(t1[m]<=78)
	{
//print(t1[m]);
	max_intensity=upper;
	}
m=m+1;
          }

         s1=derivate(derivate(t1));
          s2=derivate(derivate(t2));
          s3=derivate(derivate(t3));
          s4=derivate(derivate(t4));
          s5=derivate(derivate(t5));
          
	    



		th1=40+(FindThreshold(s1)*5);
          th2=41+(FindThreshold(s2)*5);
          th3=42+(FindThreshold(s3)*5);
          th4=43+(FindThreshold(s4)*5);
          th5=44+(FindThreshold(s5)*5);
thresh_val = newArray(5);
thresh_val[0]=th1;
thresh_val[1]=th2;
thresh_val[2]=th3;
thresh_val[3]=th4;
thresh_val[4]=th5;

count=0;
min = 220;
for (i = 0;i<5;i=i+1)
	{//print(thresh_val[i]);
if(thresh_val[i]==(40+i))
		{thresh_val[i]=0;
		count=count+1;}
else if(min>thresh_val[i]) 
		{min = thresh_val[i];
		}
	}
if (count==5){min=0;}


th = (((thresh_val[0]+thresh_val[1]+thresh_val[2]+thresh_val[3]+thresh_val[4])/5 + min)/2);
//th=min;
max_i[n-1]=max_intensity;
//print(max_i[n-1]);
if(min==0)
{threshold[n-1]=0;
//print(n+"   "+threshold[n-1]);
}
else
{threshold[n-1]=th;
avg=(avg*(denom)+th+(200-max_i[n-1]))/(denom+1);
denom++;
//print(n+"   "+avg+"   "+max_i[n-1]);
}

}			//if of iteration ends here!!!
min_diff=255;
for(n=1;n<=nSlices;n++)
	{diff=(max_i[n-1]-threshold[n-1]);
	//print("diff"+diff+"   "+threshold[n-1]);
	if(min_diff>(max_i[n-1]-threshold[n-1]))
		{min_diff=diff;}
	
	}
min_diff=max_i[default_w-1]-default_t;
for(n=1;n<=nSlices;n++)
	{if(((max_i[n-1]-threshold[n-1])>=1.5*(min_diff)||(max_i[n-1]-threshold[n-1])<=0.7*(min_diff))&&threshold[n-1]>0 )
		{t=threshold[n-1]+200-max_i[n-1];
		if(t<0){t=0;}
//print(n+"   "+avg+"   "+t);
		avg=(avg*denom-(threshold[n-1]+200-max_i[n-1]))/(denom-1);
		denom--;
		}
	}
if(denom==0){avg=0;}
avg=(avg*denom+default_t+200-max_i[default_w-1])/(denom+1);
//avg=default_t+200-max_i[default_w-1];
//print(min_diff);
suc=1;
if(avg==0){avg=155;}
//print(denom);
for (n=1; n<=nSlices; n++) 
{success=1;
if(threshold[n-1]==0||((max_i[n-1]-threshold[n-1])>=1.5*(min_diff) )||((max_i[n-1]-threshold[n-1])<=0.7*(min_diff) ))
	{
threshold[n-1]=(avg-(200-max_i[n-1]));
if(threshold[n-1]<0){threshold[n-1]=0;}
success=0;
suc=0;}
if(success==0)
{if(n>9)
{print("   "+n+"            "+floor(threshold[n-1])+"*");}
else
{print("     "+n+"            "+floor(threshold[n-1])+"*");}
}
else
{if(n>9)
{print("   "+n+"            "+floor(threshold[n-1]));}
else
{print("     "+n+"            "+floor(threshold[n-1]));}
}         
          fullname=dir+"Well"+n+"_"+name;

run("Duplicate...", "title=[I_Well"+n+"_"+name+"] duplicate range="+n+"-"+n);
run("Duplicate...", "title=[Masked_Well"+n+"_"+name+"]");

selectWindow("Masked_Well"+n+"_"+name);

setThreshold(1, threshold[n-1]);
run("Convert to Mask");
selectWindow("I_Well"+n+"_"+name);
//m=255/(max_intensity-min_intensity);
//m=255/max_intensity;
m=200-max_intensity;

//run("Subtract...", "value="+min_intensity);
//run("Multiply...", "value="+m);
//run("Invert");
//run("Add...", "value="+min_intensity);
//run("Divide...", "value=m");
//run("Subtract...", "value="+m);
run("Add...", "value="+m);
run("Invert");
saveAs("tif",dir+"I_Well"+n+"_"+name );
imageCalculator("AND", "I_Well"+n+"_"+name,"Masked_Well"+n+"_"+name);
saveAs("tif",dir+"i.m_Well"+n+"_"+name );
selectWindow("Masked_Well"+n+"_"+name);
close();
selectWindow("i.m_Well"+n+"_"+name);
//close();
//open(dir+"Well"+n+"_"+name);




 
selectWindow(name);
 }
    // rename(name);
 run("Close");
run("Images to Stack", "name=thresholded_"+name1+" title=i.m_");
run("Fire");
run("Invert LUT");
//run("Images to Stack", "method=[Copy (center)] name=thresholded_wells title=Well use");
     

  //  fullname=dir+"thresholded_"+name;




name2=replace(name, "\\.tif", "");
fullname=dir+"thresholded_"+name2+".tif";
i=0;
while(File.exists(fullname)==1)
{i=i+1;
fullname=dir+"thresholded_"+name2+"_"+i+".tif";}

//fullname=dir+"results_"+name+".txt";
//print(fullname);
saveAs("tif", fullname);







//saveAs("tif", dir+"thresholded_"+name);
  open(dir1+name1);
//selectWindow("thresholded_"+name);
selectWindow("Log");
if(suc==0){print("* : those wells where the thresholding was \n not successful using the standard method,\n threshold was scaled from other wells!");}


name1=replace(name, "\\.tif", "");
fullname=dir+"applied_thresholds_"+name1+".txt";
i=0;
while(File.exists(fullname)==1)
{i=i+1;
fullname=dir+"applied_thresholds_"+name1+"_"+i+".txt";}

//fullname=dir+"results_"+name+".txt";
//print(fullname);
saveAs("Text", fullname);



//saveAs("Text", dir+"applied_thresholds_"+name);
selectWindow("Results");
run("Close");
  
for (n=1; n<=nSlices; n++)
	{w=File.delete(dir+"I_Well"+n+"_"+name);
	 w=File.delete(dir+"i.m_Well"+n+"_"+name);
	}


}



function derivate(a) {
    b=newArray(a.length-1);
    for (n=0; n<b.length; n++) {
             b[n]=a[n+1]-a[n];
    }
    return b;
}


function PosMaxArray(a) {
    b=Array.rankPositions(a);
    return b[a.length-1];
}


function Abs(a) {
    b=newArray(a.length);
    for (n=0; nab<a.length; n++) {
           if(a[n]<0) { 
                b[n]=a[n]*(-1);
           } else {
                b[n]=a[n];
           }
    }
    return b;
}


function FindThreshold(a) {
    c=0;
t=newArray(22);
    b=PosMaxArray(a);
//max1=b;
//if(max1<=7)
//{for (x=0;x<(a.length-7);x++)
//	{
//t[x]=a[x+7];
//	}
//b=PosMaxArray(t)+7;
//max1=b;



    for (n=0; n<a.length; n++) {
          //print("b is "+a[b]);
          //print(a[b-1]);
          //print(b);
if(b<=0)
{b=0;
c=b;
n=a.length;
}

else{	
          if((a[b-1]>=a[b])||a[b]<=(0.07*a[PosMaxArray(a)])) {
               c=b;
               n=a.length;
          }
          b=b-1;
    }
    
}



return c;
}
  
      
