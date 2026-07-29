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

package ColonyArea;

import ij.*;
import ij.gui.*;
import ij.io.SaveDialog;
import ij.plugin.filter.*;
import ij.process.*;
import java.awt.*;

public class ColonyArea implements PlugInFilter {

  ImagePlus imp;

  public int setup(String arg, ImagePlus imp) {
    this.imp = imp;
    return DOES_ALL;
  }

  public int icenter(int i, double d1, double k2, double k3) {
    int c;
    c = (int) ((i - 1) * d1 + (i - 0.5) * 2 * k2 * d1 + (i - 1) * k3 * d1);
    return c;
  }

  public void run(ImageProcessor ip) {
    //(String directory = IJ.getPath(imp);
    String fname = imp.getTitle();
    IJ.showMessage(
      "Note",
      "Please select the directory and the postfix to all images."
    );
    SaveDialog sd = new SaveDialog("Save image ...", fname, "");
    String directory = sd.getDirectory();
    String fileName = sd.getFileName();
    Rectangle roi = ip.getRoi();

    ip = ip.crop();

    ImagePlus impcrop = new ImagePlus("crop", ip);

    impcrop.show();
    IJ.showMessage(
      "Note",
      "The image will be cropped and the template will be constructed."
    );

    ImageConverter iConv = new ImageConverter(impcrop);

    iConv.convertToGray8();

    impcrop.show();

    int width = ip.getWidth();
    int height = ip.getHeight();
    int w1, w2;
    int i, j;
    double del;
    double e = 5;
    double k2 = 0.06263;
    double k3 = 0.04535;
    int n1, n2;
    n1 = 0;
    n2 = 0;
    int choice = 0;
    choice =
      (int) IJ.getNumber(
        "Input \n 1 for  6 well plate. \n 2 for 12 well plate. \n 3 for 24 well plate. \n 4 for a custom plate to manually change the parameters. \n 0 to exit.",
        2
      );
    //Please enter the type of plate to be processed:\nFor 12 well plate (4X3) enter 1\nFor 24 well plate (4X6) enter 2\nFor a custom plate to change the parameters enter 3\To exit enter 0

    if (choice == 1) {
      k2 = 0.04384;
      k3 = 0.03561;
      n1 = 3;
      n2 = 2;
    } else if (choice == 2) {
      n1 = 4;
      n2 = 3;
    } else if (choice == 3) {
      k2 = 0.10387;
      k3 = 0.04277;
      n1 = 6;
      n2 = 4;
    } else if (choice == 4) {
      k2 =
        IJ.getNumber(
          "enter the ratio of well thickness to internal diameter",
          0.06263
        );
      k3 =
        IJ.getNumber(
          "enter the ratio of well spacing to internal diameter",
          0.04535
        );
      e = IJ.getNumber("enter the diameter reduction percentage", 15.0);
    } else {
      return;
    }

    w1 = (int) IJ.getNumber("enter number of wells in a row", n1);
    w2 = (int) IJ.getNumber("enter number of wells in a column", n2);

    double d1, d1w, d1h, d2, d3;

    d1w = width / (w1 + 2 * k2 * w1 + w1 * k3 - k3);
    d1h = height / (w2 + 2 * k2 * w2 + w2 * k3 - k3);
    d1 = (d1w + d1h) / 2;
    d2 = k2 * d1;
    d3 = k3 * d1;
    del = e * d1 / 100;
    double d11 = d1 - 2 * del;
    d2 = d2 + del;
    k2 = d2 / d11;
    k3 = d3 / d11;
    d1h = d1h * d11 / d1;
    d1w = d1w * d11 / d1;
    d1 = d11;
    //IJ.showMessage("ERROR!!", d1+"  "+d1h+"  "+d1w+"  "+d2+"  "+d3+"  "+del+"  "+k2+"  "+k3 );

    if ((d1w - d1h) * (d1w - d1h) > (0.02 * 0.02 * d1h * d1w)) {
      double er;
      if (d1h < d1w) {
        er = (d1w - d1h) * (w1 + 2 * w1 * k2 + w1 * k3 - k3);
        IJ.showMessage(
          "ERROR!!",
          "something is not right, please crop and straighten the image properly, or change the parameters k2 and k3. \nThe approximate max error in cropping is of the order of : " +
          er +
          "pixels horizontally."
        );
      } else {
        er = (d1h - d1w) * (w2 + 2 * w2 * k2 + w2 * k3 - k3);
        IJ.showMessage(
          "ERROR!!",
          "something is not right, please crop and straighten the image properly, or change the parameters k2 and k3. \nThe approximate max error in cropping is of the order of : " +
          er +
          "pixels vertically."
        );
      }

      return;
    }

    ImageProcessor ipt = new ByteProcessor(width, height);

    ImagePlus impt = new ImagePlus("templates_" + fileName, ipt);

    int c1, c2;
    ipt.invert();
    for (i = 1; i < w1 + 1; i++) {
      for (j = 1; j < w2 + 1; j++) {
        c1 = icenter(i, d1w, k2, k3);
        c2 = icenter(j, d1h, k2, k3);
        ipt.fillOval(c1, c2, (int) d1w, (int) d1h);
      }
    }
    ipt.invert();

    impt.show();

    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
        ip.putPixel(x, y, (ipt.getPixel(x, y) * ip.getPixel(x, y) / 255));
      }
    }
    impcrop.close();
    impt.close();
    impcrop = new ImagePlus("identified_wells_" + fileName, ip);

    iConv = new ImageConverter(impcrop);

    iConv.convertToGray8();

    impcrop.show();
    //IJ.saveAs(impcrop, "Tiff",directory+"identified_wells_"+fileName);
    IJ.wait(300);
    impcrop.close();

    ImageProcessor ipcrop = new ByteProcessor((int) d1w, (int) d1h);
    ImagePlus impwell = new ImagePlus("well", ipcrop);
    //if((w1*w2)>0){}
    ImageStack stackfinal = impwell.createEmptyStack();
    Rectangle roi1;
    //imp.close();

    for (j = 1; j < w2 + 1; j++) {
      for (i = 1; i < w1 + 1; i++) {
        c1 = icenter(i, d1w, k2, k3);
        c2 = icenter(j, d1h, k2, k3);
        ip.setRoi(c1, c2, (int) d1w, (int) d1h);
        ipcrop = ip.crop();
        stackfinal.addSlice(
          "well " + (i + w1 * (j - 1)) + " of " + w1 * w2,
          ipcrop
        );
      }
    }

    ImagePlus impwell1 = new ImagePlus("wells_" + fileName, stackfinal);
    StackConverter iConvi;
    ImageConverter iConvi_1;
    if ((w1 * w2) > 1) {
      iConvi = new StackConverter(impwell1);
      iConvi.convertToGray8();
    } //-------------------
    else {
      iConvi_1 = new ImageConverter(impwell1);
      iConvi_1.convertToGray8();
    }

    impwell1.show();

    IJ.saveAs(impwell1, "Tiff", directory + "wells_" + fileName);

    IJ.run(impwell1, "Colony thresolder", "");
  }
}
