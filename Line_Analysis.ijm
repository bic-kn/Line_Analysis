// BIC_Recruitment in line ROI
// =======================================================
// Copyright 2015 BioImaging Center, University of Konstanz
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version. 
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. If data produced 
// with this program or a derivative of this program is used for 
// publications the original authors should be acknowledged appropriately. 

// See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// =======================================================

// BioImaging Center, University of Konstanz <bioimaging@uni-konstanz.de>
//   Martin St�ckl <martin.stoeckl@uni-konstanz.de>

// =======================================================

	
	var VOID = 0;
	var row = 0;	
	var ROInum = 0;
	var chanum = 0;
	var endimage = 0;
	var rtime = 0;
	var widthline =0;
	var heightline = 0;
 	var filename = "";
	var path ="";
	offset = 14;
	
macro "Line analysis staining" {
	
	method = newArray("line ROI", "whole nucleus", "single slide");

	run("Close All");	
	run("Clear Results");
	roiManager("Reset");
	
	Dialog.create("Method selection");
		Dialog.addCheckbox("Run in batch mode?", true);
		Dialog.addChoice("Select bgd selection method:", method, method[1]);
		Dialog.addCheckbox("Multichannel measurement?", false);
		Dialog.show();
	batch = Dialog.getCheckbox(); 
	method[0] = Dialog.getChoice();
	multi_ch = Dialog.getCheckbox(); 

	if (multi_ch && method[0] != "single slide")
	
		exit("At the moment the macro supports multichannels only for single slide measurements!");

		
	file = File.openDialog("Select first image");
	if (endsWith(file, ".lsm")) {
		run("Show LSMToolbox","ext");
		Ext.lsmOpen(file);
	}
	else
	 {
		run("Bio-Formats (Windowless)", "open=["+file+"]");
		if (method[0] != "single slide") 
			exit("At the moment the macro supports only .lsm/nfor time series measurements!");
	}	
	files=getTitle();
//	prefix=substring(files, 0, lengthOf(files)-offset);
	path = getDirectory("image") + File.separator;
	
	if (batch) {
		Dialog.create("Select the prefix");
			Dialog.addString("Prefix used for image selection:", files, lengthOf(files));
		Dialog.show();
			prefix = Dialog.getString();

		list = getFileList(path);
		files = newArray();
		for (i=0; i<list.length; i++) 
//			if (startsWith(list[i], prefix) && endsWith(list[i], ".lsm")) 
			if (startsWith(list[i], prefix)) 
				files = Array.concat(files, list [i]);
		Array.sort(files);

		Dialog.create("Image files detected!");
			Dialog.addMessage(""+files.length+"  files will be analyzed!");
			Dialog.show();
	}
	
	
	if (multi_ch) {
		
		excit = newArray("405 nm", "458 nm", "488 nm", "514 nm", "543 nm", "633 nm", "TL", "Select");
		getDimensions(VOID, VOID, channels, VOID, VOID);
		choice = newArray(channels);	
	
		Dialog.create("Wavelength selection");
			for (i=0; i<channels;i++) 
				Dialog.addChoice("Excitation wavelength of Channel "+ i+1 +":", excit, "Select");
				Dialog.addChoice("Do thresholding in Channel:", excit, "633 nm");
		Dialog.show;
			for (i=0; i<channels;i++) 	
				choice[i] = Dialog.getChoice();
			chan_thresh = Dialog.getChoice();

		if (chan_thresh == "Select")
			chan_thresh = 1;
		else {
			i = 0;
			while (chan_thresh != choice[i] && i < choice.length)
				i++;
			chan_thresh = i+1;
		}		
		chan_sel = Channel_Select(choice);
		chan_nums =newArray();
		for (i=1; i<=choice.length;i++) 
			if (choice[i-1] != "Select") 
				chan_nums = Array.concat(chan_nums, i);
		
		for (i=0; i<choice.length;i++) 
			if (choice[i] == "Select")
				choice = Array_Remove(choice, i);
		
	}
	
	if (File.exists(path + "line.txt")) {
		linebox = File.openAsString(path+"line.txt");
		linebox = split(linebox, ",");
		for (i=0; i<linebox.length; i++) 
			linebox[i] = parseInt(linebox[i]);
		widthline = linebox[0];
		heightline = linebox[1];
		chanum = linebox[2];
		xsel = linebox[3];
		ysel = linebox[4];
		makeRectangle(xsel, ysel, widthline, heightline);
		Stack.setChannel(chanum);
		waitForUser("Check line size and position! To change for the series delete line.txt");
	}
	else {
		waitForUser("Select the channel to measure and mark line with box then press OK!");
			getSelectionBounds(xsel, ysel, widthline, heightline);
			Stack.getPosition(chanum, VOID, VOID);
		File.saveString(widthline+","+heightline+","+chanum+","+xsel+","+ysel, path+"line.txt");
	}

	
	if (multi_ch)
 {
		chanum_old = chanum;
		for(i=0; i<chanum_old-1;i++) 
			if(chan_sel[i] == 0)
				chanum--;
				
		thres_old = chan_thresh;
		for(i=0; i<thres_old-1;i++) 
			if(chan_sel[i] == 0)
				chan_thresh--;

	}
	else
		chan_sel = toString(chanum);
			
	if (method[0] != "single slide") {
		pb = 1;
		irr = 1;
		rtime = 60;
		Stack.getUnits(VOID, VOID, VOID, Time, VOID);
		getDimensions(VOID, VOID, VOID, VOID, frames);	
		
		do {
			 
			Dialog.create("Parameter setup");
				Dialog.addNumber("Number of not irradiated images", pb);
				Dialog.addNumber("Number of irradiation images", irr);
				Dialog.addNumber("Recruitment time after irradiation [in "+Time+"]", rtime);	
				Dialog.show();
			pb = Dialog.getNumber();
			irr = Dialog.getNumber();
			solltime = Dialog.getNumber();
			
			tstamp = Tstamps(file);
						
			if (tstamp[frames-1] - solltime < tstamp[pb+irr-1]) {
				rtime = tstamp[frames-1] - tstamp[pb+irr-1];
				Dialog.create("ERROR");
					Dialog.addCheckbox("Time series is too short for selected recruitment time!\n Continue anyway?", false);
					Dialog.addMessage("Approx. length of time series: "+rtime+" "+Time);
					Dialog.show();
				control = Dialog.getCheckbox(); 
				solltime = rtime;
			}
			else 		
				control = 1;
						
		} while (control == 0);
	}
	File.makeDirectory(path+"Results");
	
	
	if (batch == 1) {
		run("Close All");
		rtime = newArray(files.length);
		endimage = newArray(files.length);

		if (method[0] == "whole nucleus")
			for (i=files.length; i>0; i--) {
				filename = path+files[i-1];
				Ext.lsmOpen(filename);
				tstamp = Tstamps(filename);
				endimage[i-1] = findTime(solltime, tstamp, irr);
				rtime[i-1] = tstamp[endimage[i-1]-1] - tstamp[pb+irr-1];
				ProcessNuc(files[i-1], endimage[i-1]);
			}

		if (method[0] == "line ROI")
			for (i=files.length; i>0; i--) {
				filename = path+files[i-1];
				Ext.lsmOpen(filename);
				tstamp = Tstamps(filename);
				endimage[i-1] = findTime(solltime, tstamp, irr);
				rtime[i-1] = tstamp[endimage[i-1]-1] - tstamp[pb+irr-1];
				ProcessLine(files[i-1], endimage[i-1]);
			}

		if (method[0] == "single slide")
			for (i=files.length; i>0; i--) {
				filename = path+files[i-1];
				run("Bio-Formats (Windowless)", "open=["+filename+"]");
				Stack.setChannel(chanum_old);
				ProcessNuc_single(files[i-1], chan_sel, chanum, chan_thresh);
			}

		run("Clear Results");
		if (method[0] != "single slide")
 		{
			for (i=files.length; i>0; i--)
				Analyse(files[i-1], rtime[i-1], endimage[i-1]);
			saveAs("Results", path+"Results"+ File.separator + prefix+"_Results.txt");
		}
		else
 		{
			for (j=0; j<choice.length; j++) {
				ROInum = 0;
				row = 0;
				run("Clear Results");
				for (i=files.length; i>0; i--)
 
					Analyse_single(files[i-1], j+1);
				saveAs("Results", path+"Results"+ File.separator + prefix+"_"+choice[j]+"_Results.txt");
			}
		}
	}
	else {
		if (method[0] != "single slide") {
			endimage = findTime(solltime, tstamp, irr);
			rtime = tstamp[endimage-1] - tstamp[pb+irr-1];
		
			if (method[0] == "whole nucleus")
				ProcessNuc(files);
			else 
				ProcessLine(files);
			run("Clear Results");
			Analyse(files, rtime);
			saveAs("Results", path+"Results"+ File.separator + files +"_Results.txt");
		}
		else {
			ProcessNuc_single(files, chan_sel, chanum, chan_thresh);
			for (j=0; j<chan_nums.length; j++) {
				run("Clear Results");
				ROInum = 0;
				row = 0;
				Analyse_single(files, j+1);
				saveAs("Results", path+"Results"+ File.separator + files +"_"+choice[j]+"_Results.txt");
			}
		}
			
	}

	if (batch)
		roiManager("Save", path+"Results"+ File.separator + prefix+"_ROISet.zip");
	else
		roiManager("Save", path+"Results"+ File.separator + files+"_ROISet.zip");

}


	
function Tstamps(fullname) {
	stamp = split(Ext.getTStamps(fullname), ",");
	for (i=0; i<lengthOf(stamp); i++) 
		stamp[i] = parseFloat(stamp[i]);
	return stamp;
}

function findTime(time, stamp, laserimg) {

	tempstamp = newArray(stamp.length);
			for (i = 0; i<stamp.length; i++) {
				tempstamp[i] = stamp[i] - stamp[pb+laserimg-1] - time;
				if (tempstamp[i] < 0)
					tempstamp[i] *= -1;
			 }
			Array.getStatistics(tempstamp, min, VOID, VOID, VOID);
		i = 0;
		while (tempstamp[i] - min != 0) 
			i++;
		image = i + 1;
		return image;
		
}

function ProcessNuc(name, end_image)  {
	run("Enhance Contrast", "saturated=0.35");
	waitForUser("Mark cell with box and press OK!");
	getSelectionBounds(box_x, box_y, VOID, VOID);
	run("Duplicate...", "title="+name+"_crop duplicate channels="+chanum+" frames=1-"+frames);
	run("Duplicate...", "title=blur duplicate channels="+chanum+" frames=1-"+frames);
	selectWindow(name);
 	close();
	selectWindow(name+"_crop");
	save(path+"Results"+ File.separator + name+".tiff");
	getDimensions(imgwidth, imgheight, VOID, VOID, VOID);
	Stack.setFrame(end_image);
	xsel_temp = xsel - box_x;
	ysel_temp = ysel - box_y;
	if (xsel_temp > 0 && ysel_temp > 0) {
		getDimensions(check_w, check_h, VOID, VOID, VOID);
			if (xsel_temp + widthline >= check_w || ysel_temp + heightline >= check_h)
				makeRectangle(2, 2, widthline, heightline);
			else
				makeRectangle(xsel_temp, ysel_temp, widthline, heightline);
	}
	else
		makeRectangle(2, 2, widthline, heightline);

	waitForUser("Move box on line!");
	
	run("Add to Manager");
	Overlay.hide;
	roiManager("Select", roiManager("count")-1);
	roiManager("Rename", "line_"+name);
	ROIbefore = roiManager("count");
	selectWindow("blur");
	run("Gaussian Blur...", "sigma=3 stack");
	run("Threshold...");
	setAutoThreshold("Huang dark");
	Stack.setChannel(chanum);
	Stack.setFrame(end_image);
	
	waitForUser("Adjust threshold and press OK!");
	run("Analyze Particles...", "size=100-Infinity circularity=0.00-1.00 show=Nothing exclude include add slice");	
	ROIafter = roiManager("count");
	while (ROIafter-ROIbefore > 1) {
		waitForUser("More than 1 ROI created! Remove unwanted ROI(s) and press OK!");
		ROIafter = roiManager("count");
	}
	while (ROIafter-ROIbefore == 0) {
		setTool("freehand");
		waitForUser("No ROI created! Create selection by hand and press OK!");
		run("Add to Manager");
		ROIafter = roiManager("count");
		setTool("rectangle");
		}
	ROIafter = roiManager("count");
	roiManager("Select", roiManager("count")-1);
	roiManager("Rename", "nuc_"+name);
	Overlay.hide;
	close();
}

function ProcessNuc_single(name, chan_selection, channel_number, thresh_number)  {
	run("Enhance Contrast", "saturated=0.35");
	waitForUser("Mark cell with box and press OK!");
	getSelectionBounds(box_x, box_y, VOID, VOID);
	run("Duplicate...", "title="+name+"_crop duplicate channels=1-"+channels);
	Slice_Remove(chan_selection);
	getDimensions(VOID, VOID, channels, VOID, VOID);
	run("Duplicate...", "title=blur duplicate channels=1-"+channels);
	selectWindow(name);
 	close();
	selectWindow(name+"_crop");
	save(path+"Results"+ File.separator + name+".tiff");
	getDimensions(imgwidth, imgheight, VOID, VOID, VOID);
	xsel_temp = xsel - box_x;
	ysel_temp = ysel - box_y;
	if (xsel_temp > 0 && ysel_temp > 0) {
		getDimensions(check_w, check_h, VOID, VOID, VOID);
			if (xsel_temp + widthline >= check_w || ysel_temp + heightline >= check_h)
				makeRectangle(2, 2, widthline, heightline);
			else
				makeRectangle(xsel_temp, ysel_temp, widthline, heightline);
	}
	else
		makeRectangle(2, 2, widthline, heightline);

	Stack.setChannel(channel_number);
	waitForUser("Move box on line!");
	run("Add to Manager");
	roiManager("Select", roiManager("count")-1);
	roiManager("Rename", "line_"+name);
	Overlay.hide;
	waitForUser("Select background!");
	run("Add to Manager");
	roiManager("Select", roiManager("count")-1);
	roiManager("Rename", "bgd_"+name);
	Overlay.hide;
	ROIbefore = roiManager("count");
	selectWindow("blur");
	run("Gaussian Blur...", "sigma=3 stack");
	Stack.setChannel(thresh_number);
	setAutoThreshold("Huang dark");
	run("Threshold...");
	
	waitForUser("Adjust threshold and press OK!");
	run("Analyze Particles...", "size=100-Infinity circularity=0.00-1.00 show=Nothing exclude include add slice");	
	ROIafter = roiManager("count");
	while (ROIafter-ROIbefore > 1) {
		waitForUser("More than 1 ROI created! Remove unwanted ROI(s) and press OK!");
		ROIafter = roiManager("count");
	}
	while (ROIafter-ROIbefore == 0) {
		setTool("freehand");
		waitForUser("No ROI created! Create selection by hand and press OK!");
		run("Add to Manager");
		ROIafter = roiManager("count");
		setTool("rectangle");
		}
	ROIafter = roiManager("count");
	roiManager("Select", roiManager("count")-1);
	roiManager("Rename", "nuc_"+name);
	somanyROIs = roiManager("count")-3;
	roiManager("Select", newArray(somanyROIs,somanyROIs+2));
	roiManager("XOR");
	roiManager("Add");
	roiManager("Select", roiManager("count")-1);
	roiManager("Rename", "nucOnly_"+name);
	Overlay.hide;
	close();
}



function ProcessLine(name, end_image)  {
	run("Enhance Contrast", "saturated=0.35");
	waitForUser("Mark cell with box and press OK!");
	run("Duplicate...", "title="+name+"_crop duplicate channels="+chanum+" frames=1-"+frames);
	selectWindow(name);
 	close();
	selectWindow(name+"_crop");
	save(path+"Results"+ File.separator + name+".tiff");
	makeRectangle(2, 2, widthline, heightline);
	Stack.setFrame(end_image);
	waitForUser("Move box on line!");
	run("Add to Manager");
	roiManager("Select", roiManager("count")-1);
	roiManager("Rename", "line_"+name);
	getSelectionBounds(xpos, ypos, widthline, heightline);
	makeRectangle(xpos+widthline+5, ypos, widthline, heightline);
	waitForUser("Move box to nucleus!");	
	run("Add to Manager");
	roiManager("Select", roiManager("count")-1);
	roiManager("Rename", "nuc_bgd_"+name);
	Overlay.hide;
}
				
function Analyse(name, time, end_image) {
	selectWindow(name+"_crop");
	Stack.setFrame(end_image);
	roiManager("Select", ROInum);
	Stack.setFrame(end_image);
	getStatistics(linearea, linemean, VOID, max);
	cont = 1;
	if (max == 4095)													// imported .lsm are treated as 16-bit images regardless of real bit depth, adjust
		cont = getBoolean("Overexposed pixels detected! Use image anyway?");
	if (cont) {
		getSelectionBounds(line_x, line_y, VOID, VOID);
		setResult("Label", row, name);
		setResult("Slide #", row, end_image);
		setResult("Recruit. Time", row, time);
		setResult("Line area", row, linearea);	
		setResult("Line mean ", row, linemean);
		ROInum++;
		roiManager("Select", ROInum);
		Stack.setFrame(end_image);
		getStatistics(nucarea, nucmean);
		setResult("Nucleus area", row, nucarea);
		setResult("Nucleus mean ", row, nucmean);
		getSelectionBounds(nuc_x, nuc_y, VOID, VOID);
		Overlay.hide;
		waitForUser("Select background!");
		run("Add to Manager");
		roiManager("Select", roiManager("count")-1);
		roiManager("Rename", "bgd_"+name);		
		getStatistics(VOID, bgd);
		roiManager("Select", ROInum);
		Stack.setFrame(pb);
		waitForUser("Adjust nucleus location!");	
		run("Add to Manager");
		roiManager("Select", roiManager("count")-1);
		roiManager("Rename", "nuc_preirr_"+name);
		getSelectionBounds(prenuc_x, prenuc_y, VOID, VOID);
		getStatistics(VOID, prenuc);
		setResult("Nuc. Preirr. mean ", row, prenuc);
		makeRectangle(prenuc_x - nuc_x + line_x, prenuc_y - nuc_y + line_y, widthline, heightline);
		run("Add to Manager");
		roiManager("Select", roiManager("count")-1);
		roiManager("Rename", "line_preirr_"+name);
		getStatistics(VOID, preline);
		setResult("Line Preirr. mean", row, preline);
		setResult("Background", row, bgd);
		percent_nuc = (nucmean-bgd) / (prenuc-bgd) * 100;
		subtr_line = linemean  * (prenuc-bgd) / (nucmean-bgd) - preline;
		corr_line = 100 * ((linemean-bgd) * (prenuc-bgd)) / ((preline-bgd) * (nucmean-bgd));
		setResult("Nucleus [Percent of preirr. nucleus]", row, percent_nuc);
		setResult("Line [Percent of preirr. line]", row, corr_line);
		setResult("Line [Int. - preirr. int.]", row, subtr_line);
		updateResults();
		row++;
		ROInum++;
	}
	else
		ROInum += 2;
}

function Analyse_single(name, channel_number) {

		selectWindow(name+"_crop");
		roiManager("Select", ROInum);
		Stack.setChannel(channel_number);
		getStatistics(linearea, linemean, VOID, max);
		setResult("Label", row, name);
		setResult("Line area", row, linearea);	
		setResult("Line mean ", row, linemean);
		ROInum++;
		roiManager("Select", ROInum);
		Stack.setChannel(channel_number);
		getStatistics(VOID, bgd);
		setResult("Background", row, bgd);
		ROInum++;
		roiManager("Select", ROInum);
		Stack.setChannel(channel_number);
		getStatistics(nucarea, nucmean);
		setResult("Nucleus area", row, nucarea);
		setResult("Nucleus mean ", row, nucmean);
		ROInum++;
		roiManager("Select", ROInum);
		Stack.setChannel(channel_number);
		getStatistics(nucOnlyarea, nucOnlymean,VOID,VOID,nucOnlydev);
		setResult("Nucleus area w/o line", row, nucOnlyarea);
		setResult("Nucleus mean w/o line", row, nucOnlymean);
		setResult("Nucleus  w/o line STDDEV", row, nucOnlydev);
		nuc_corr = nucmean - bgd;
		line_corr = linemean  - bgd;
		nucOnly_corr = nucOnlymean - bgd;
		setResult("Nucleus [Nucleus mean corr.]", row, nuc_corr);
		setResult("Nucleus w/o line [Nucleus mean corr.]", row, nucOnly_corr);		
		setResult("Line [Line mean corr.]", row, line_corr);
		//setResult("Total intensity [whole nucleus]", row, nuc_corr*nucarea);
		updateResults();
		
		row++;
		ROInum++;
}	
/*		 
function Channel_Select (sel) {

	chan ="";
	for (i=1; i<=sel.length;i++) 
		if (sel[i-1] != "Select") 
			chan+=toString(i)+",";
	chan = substring(chan, 0, lengthOf(chan)-1);
	return chan;
}

*/

function Array_Remove(a,pos) {
	a1 = Array.slice(a, 0, pos);
	a2 = Array.slice(a, pos+1);
	a = Array.concat(a1, a2);
	return a;
}
								

function Channel_Select(chan) {
	sel = newArray();
	for (i=0; i<chan.length; i++) {
		if (chan[i] == "Select") 
			sel = Array.concat(sel, 0);
		else
			sel = Array.concat(sel, 1);
	}
	return sel;
}
	
function Slice_Remove(chan) {
	for (i=channels; i>0; i--) {
		Stack.setChannel(i);
		if (chan[i-1] == 0)
			run("Delete Slice", "delete=channel");
	}
}

