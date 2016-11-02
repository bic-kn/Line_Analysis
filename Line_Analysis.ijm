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
// Martin Stoeckl <martin.stoeckl@uni-konstanz.de>
// Stefan Helfrich <stefan.helfrich@uni-konstanz.de>
// =======================================================

// TODO Remove global variables where possible

var VOID = 0;
var row = 0;	
var ROInum = 0;
var chanum = 0;
var endimage = 0;
var recruitmentTime = 0;
var widthline = 0;
var heightline = 0;
var filename = "";
var path ="";
offset = 14;

macro "Line analysis staining" {

	
	// Reset everything
	run("Close All");	
	run("Clear Results");
	roiManager("Reset");
	
	// Create and show dialog
	Dialog.create("Method selection");
	Dialog.addCheckbox("Batch mode", true);
	methods = newArray("line ROI", "whole nucleus", "single slide");
	Dialog.addChoice("Background selection method:", methods, methods[1]);
	Dialog.addCheckbox("Multichannel measurement", false);
	Dialog.show();

	// Get user input
	batch = Dialog.getCheckbox(); 
	method = Dialog.getChoice();
	multi_ch = Dialog.getCheckbox(); 

	if (multi_ch && method != "single slide")
 {
		exit("Multichannels only supported for single slide measurements");
	}

	file = File.openDialog("Select first image");
	if (endsWith(file, ".lsm")) {
		// Use LSM Toolbox for opening LSM instead of Bio-Formats
		run("Show LSMToolbox", "ext");
		Ext.lsmOpen(file);
	}
 else {
		// Use Bio-Formats for other formats
		run("Bio-Formats (Windowless)", "open=["+file+"]");
		if (method != "single slide") 
{
			exit("Only .lsm/nfor for time series measurements supported");
		}
	}	
	files = getTitle();
//	prefix=substring(files, 0, lengthOf(files)-offset);
	path = getDirectory("image") + File.separator;

	if (batch) {
		// Dialog for prefix selection
		Dialog.create("Select the prefix");
		Dialog.addString("Prefix used for image selection:", files, lengthOf(files));
		Dialog.show();
		prefix = Dialog.getString();

		// Filter file list according to prefix
		list = getFileList(path);
		files = newArray();
		for (i=0; i<list.length; ++i) 
{
//			if (startsWith(list[i], prefix) && endsWith(list[i], ".lsm")) 
			if (startsWith(list[i], prefix)) { 
				files = Array.concat(files, list[i]);
			}
		}
		Array.sort(files);

		showMessage("Image files detected!", "&files.length files will be analyzed!");
	}
	
	// Get user input on multi-channel measurements
	if (multi_ch) {
		excit = newArray("405 nm", "458 nm", "488 nm", "514 nm", "543 nm", "633 nm", "TL", "Select");

		// Get number of channels
		getDimensions(VOID, VOID, channels, VOID, VOID);
	
		// Create dialog
		Dialog.create("Wavelength selection");
		for (i=0; i<channels; ++i) { 
			Dialog.addChoice("Excitation wavelength of Channel "+ i+1 +":", excit, "Select");
		}
		Dialog.addChoice("Do thresholding in Channel:", excit, "633 nm");
		Dialog.show();

		// Get input from dialog
		excitationWavelengths = newArray(channels);
		for (i=0; i<channels; ++i) {
			excitationWavelengths[i] = Dialog.getChoice();
		}
		chan_thresh = Dialog.getChoice();
 // String representation

		if (chan_thresh == "Select")
 {
			chan_thresh = 1; // Default channel
		} else {
			chan_thresh = convertWavelengthToChannelNumber(chan_thresh, exicationWavelengths);
		}

		// Create an array of channels for which wavelengths have been selected
		chan_sel = Channel_Select(excitationWavelengths);
		chan_nums = newArray();
		for (i=0; i<excitationWavelengths.length; ++i) { 
			if (excitationWavelengths[i] != "Select") 
{
				chan_nums = Array.concat(chan_nums, i+1);
			}
		}

		// Overwrite excitationWavelenghts s.t. it only contains specified wavelengths 
		for (i=0; i<excitationWavelengths.length; ++i) { 
			if (excitationWavelengths[i] == "Select")
 {
				excitationWavelengths = Array_Remove(excitationWavelengths, i);
			}
		}
	}
	
	// Read line for evaluation from file
	if (File.exists(path + "line.txt")) {
		// Load selection from file
		makeSelectionFromFile(path+"line.txt");

		// Let the user validate the selection
		waitForUser("Check line size and position! To change for the series delete line.txt");
	}
 else {
		// Let the user select a rectangle
		waitForUser("Select the channel to measure and mark line with box then press OK!");

		// Save selection to file
		saveSelectionToFile(path+"line.txt");
	}

	if (multi_ch)
 {
		// Update chanum
		Stack.getPosition(chanum, VOID, VOID);
		chanum_old = chanum;
		for(i=0; i<chanum_old-1; ++i) 
{
			if(!chan_sel[i])
 {
				chanum--;
			}
		}

		// Update chan_thresh
		thres_old = chan_thresh;
		for (i=0; i<thres_old-1; ++i) {
			if (!chan_sel[i]) {
				chan_thresh--;
			}
		}
	}
 else
 {
		chan_sel = toString(chanum);
	}
			
	if (method != "single slide") {
		// Number of frames before bleaching
		beforeBleachingFrames = 1;
		// Number of frames that are not analysed because the images are basically white
		irradiationFrames = 1;
		recruitmentTime = 60;
		Stack.getUnits(VOID, VOID, VOID, timeUnit, VOID);
		getDimensions(VOID, VOID, VOID, VOID, frames);

		// Repeatedly ask the user for input and verify the input. If the inputs are good, continue.
		do {
			Dialog.create("Parameter setup");
			Dialog.addNumber("Number of not irradiated images", beforeBleachingFrames);
			Dialog.addNumber("Number of irradiation images", irradiationFrames);
			Dialog.addNumber("Recruitment time after irradiation [in &timeUnit]", recruitmentTime);	
			Dialog.show();
			
			beforeBleachingFrames = Dialog.getNumber();
			irradiationFrames = Dialog.getNumber();
			analysisTime = Dialog.getNumber();

			// Get timestamps from file
			timeStamps = Tstamps(file);
			
			// Check if experiment is long enough to cover the recruitment
			if (timeStamps[frames-1] - timeStamps[beforeBleachingFrames+irradiationFrames-1] < analysisTime) {
				recruitmentTime = timeStamps[frames-1] - timeStamps[beforeBleachingFrames+irradiationFrames-1];
				
				Dialog.create("ERROR");
				Dialog.addCheckbox("Time series is too short for selected recruitment time!\n Continue anyway?", false);
				Dialog.addMessage("Approximated length of time series: &recruitmentTime &timeUnit");
				Dialog.show();
				
				control = Dialog.getCheckbox(); 
				analysisTime = recruitmentTime;
			}
 else { 		
				control = true;
			}			
		} while (!control);
	}
	
	// Create results directory
	File.makeDirectory(path+"Results");
	
	if (batch) {
		/*
		 * Batch mode
		 */
		run("Close All");
		recruitmentTime = newArray(files.length);
		endimage = newArray(files.length);

		// TODO Refactor similar code
		if (method == "whole nucleus")
 {
			for (i=files.length; i>0; i--) {
				filename = path+files[i-1];
				Ext.lsmOpen(filename);
				timeStamps = Tstamps(filename);
				endimage[i-1] = findTime(analysisTime, timeStamps, irradiationFrames);
				recruitmentTime[i-1] = timeStamps[endimage[i-1]-1] - timeStamps[beforeBleachingFrames+irradiationFrames-1];
				
				ProcessNuc(files[i-1], endimage[i-1]);
			}
		}

		if (method == "line ROI")
 {
			for (i=files.length; i>0; i--) {
				filename = path+files[i-1];
				Ext.lsmOpen(filename);
				timeStamps = Tstamps(filename);
				endimage[i-1] = findTime(analysisTime, timeStamps, irradiationFrames);
				recruitmentTime[i-1] = timeStamps[endimage[i-1]-1] - timeStamps[beforeBleachingFrames+irradiationFrames-1];
				
				ProcessLine(files[i-1], endimage[i-1]);
			}
		}

		if (method == "single slide")
 {
			for (i=files.length; i>0; i--) {
				filename = path+files[i-1];
				run("Bio-Formats (Windowless)", "open=["+filename+"]");
				Stack.setChannel(chanum_old);
				
				ProcessNuc_single(files[i-1], chan_sel, chanum, chan_thresh);
			}
		}

		run("Clear Results");

		if (method != "single slide")
 {
			for (i=files.length; i>0; i--)
 {
				Analyse(files[i-1], recruitmentTime[i-1], endimage[i-1]);
			}
			saveAs("Results", path+"Results"+File.separator+prefix+"_Results.txt");
		}
 else {
			for (j=0; j<choice.length; j++) {
				ROInum = 0;
				row = 0;
				run("Clear Results");
				for (i=files.length; i>0; i--)
 {
 
					Analyse_single(files[i-1], j+1);
				}
				
				saveAs("Results", path+"Results"+ File.separator + prefix+"_"+choice[j]+"_Results.txt");
			}
		}
	}
 else {
		/* 
		 * Non-batch mode
		 */
		if (method != "single slide") {
			// TODO ???
			endimage = findTime(analysisTime, timeStamps, irradiationFrames);
			recruitmentTime = timeStamps[endimage-1] - timeStamps[beforeBleachingFrames+irradiationFrames-1];

			if (method == "whole nucleus")
 {
				ProcessNuc(files);
			} else { 
				ProcessLine(files);
			}

			run("Clear Results");

			Analyse(files, recruitmentTime);
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
 {
		roiManager("Save", path+"Results"+ File.separator + prefix+"_ROISet.zip");
	} else {
		roiManager("Save", path+"Results"+ File.separator + files+"_ROISet.zip");
	}
}


/*
 * Get timestamps from metadata.
 */
function Tstamps(fullname) {
	stamp = split(Ext.getTStamps(fullname), ",");
	for (i=0; i<lengthOf(stamp); i++) 
{
		stamp[i] = parseFloat(stamp[i]);
	}
	return stamp;
}

/*
 * Find the frame number of the last frame s.t. the analysis time is best covered.
 */
function findTime(analysisTime, timeStamps, laserimg) {
	// TODO pb???
	tempstamp = newArray(timeStamps.length);
	analysisStartTime = timeStamps[pb+laserimg-1];
	for (i = 0; i<timeStamps.length; i++) {
		tempstamp[i] = abs((timeStamps[i] - analysisStartTime) - analysisTime);
	}

	// Get the index of the minimum value
	ranks = Array.rankPositions(tempstamp);

	return ranks[0] + 1;
		
}

/*
 * TODO Documentation
 */
function ProcessNuc(name, end_image)  {
	run("Enhance Contrast", "saturated=0.35");
	
	waitForUser("Mark cell with box and press OK!");

	// Crop active image
	getSelectionBounds(box_x, box_y, VOID, VOID);
	run("Duplicate...", "title="+name+"_crop duplicate channels="+chanum+" frames=1-"+frames);
	
	// Create duplicate for blurring
	run("Duplicate...", "title=blur duplicate channels="+chanum+" frames=1-"+frames);

	// Close original image
	selectWindow(name);
 	close();

 	// Save cropped image to results folder
	selectWindow(name+"_crop");
	save(path+"Results"+File.separator+name+".tiff");

	// Shift ROI in cropped image
	getDimensions(imgwidth, imgheight, VOID, VOID, VOID);
	Stack.setFrame(end_image);
	xsel_temp = xsel - box_x;
	ysel_temp = ysel - box_y;
	if (xsel_temp > 0 && ysel_temp > 0) {
		getDimensions(check_w, check_h, VOID, VOID, VOID);
			if (xsel_temp + widthline >= check_w || ysel_temp + heightline >= check_h)
 {
				makeRectangle(2, 2, widthline, heightline);
			} else {
				makeRectangle(xsel_temp, ysel_temp, widthline, heightline);
			}
	}
 else {
		makeRectangle(2, 2, widthline, heightline);
	}

	// Ask user for validation
	waitForUser("Move box on line!");

	// Add to ROI manager
	run("Add to Manager");
	
	Overlay.hide();
	roiManager("Select", roiManager("count")-1);
	roiManager("Rename", "line_"+name);
	ROIbefore = roiManager("count");
	
	// Blur image
	selectWindow("blur");
	run("Gaussian Blur...", "sigma=3 stack");
	
	// Threshold image
	run("Threshold...");
	setAutoThreshold("Huang dark");
 // default threshold
	Stack.setChannel(chanum);
	Stack.setFrame(end_image);

	// Ask user for validation
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

/*
 * TODO Documentation
 */
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


/*
 * TODO Documentation
 */
function ProcessLine(name, end_image) {
	// Enhance contrast
	run("Enhance Contrast", "saturated=0.35");
	
	// Let the user select the cell
	waitForUser("Mark cell with box and press OK!");

	// Crop active image
	run("Duplicate...", "title="+name+"_crop duplicate channels="+chanum+" frames=1-"+frames);
	
	// Close original image
	selectWindow(name);
 	close();

 	// Save cropped image to results folder
	selectWindow(name+"_crop");
	save(path+"Results"+ File.separator + name+".tiff");

	// Propose line selection
	makeRectangle(2, 2, widthline, heightline);
	Stack.setFrame(end_image);
	waitForUser("Move box on line!");

	// Add line selection to ROI Manager
	run("Add to Manager");
	roiManager("Select", roiManager("count")-1);
	roiManager("Rename", "line_"+name);
	
	// Propose selection for nucleus
	getSelectionBounds(xpos, ypos, widthline, heightline);
	makeRectangle(xpos+widthline+5, ypos, widthline, heightline);
	waitForUser("Move box to nucleus!");
	
	// Add nucleus selection
 to ROI Manager
	run("Add to Manager");
	roiManager("Select", roiManager("count")-1);
	roiManager("Rename", "nuc_bgd_"+name);
	Overlay.hide;
}

/*
 * TODO Documentation
 */
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
		// Write information about line to results table
		getSelectionBounds(line_x, line_y, VOID, VOID);
		setResult("Label", row, name);
		setResult("Slide #", row, end_image);
		setResult("Recruit. Time", row, time);
		setResult("Line area", row, linearea);	
		setResult("Line mean ", row, linemean);
		ROInum++;

		// Write information about nucleus to results table
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
 else {
		ROInum += 2;
	}
}

/*
 * TODO Documentation
 * TODO Remove ROInum (global variable)
 * TODO Remove row (global variable)
 */
function Analyse_single(name, channel_number) {
		selectWindow(name+"_crop");

		// Get info about line
		roiManager("Select", ROInum);
		Stack.setChannel(channel_number);
		getStatistics(linearea, linemean, VOID, max);
		setResult("Label", row, name);
		setResult("Line area", row, linearea);	
		setResult("Line mean ", row, linemean);
		ROInum++;

		// Get info about background
		roiManager("Select", ROInum);
		Stack.setChannel(channel_number);
		getStatistics(VOID, bgd);
		setResult("Background", row, bgd);
		ROInum++;

		// Get info about nucleus
		roiManager("Select", ROInum);
		Stack.setChannel(channel_number);
		getStatistics(nucarea, nucmean);
		setResult("Nucleus area", row, nucarea);
		setResult("Nucleus mean ", row, nucmean);
		ROInum++;

		// Get info about nucleus without line
		roiManager("Select", ROInum);
		Stack.setChannel(channel_number);
		getStatistics(nucOnlyarea, nucOnlymean,VOID,VOID,nucOnlydev);
		setResult("Nucleus area w/o line", row, nucOnlyarea);
		setResult("Nucleus mean w/o line", row, nucOnlymean);
		setResult("Nucleus  w/o line STDDEV", row, nucOnlydev);

		// Compute results
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
 * Remove element from array.
 */
function Array_Remove(a, pos) {
	a1 = Array.slice(a, 0, pos);
	a2 = Array.slice(a, pos+1);
	a = Array.concat(a1, a2);
	return a;
}
								

/*
 * Converts from a string array to a boolean area by checking the elements for equality to "Select"
 */
function Channel_Select(chan) {
	channelSelected = newArray();
	for (i=0; i<chan.length; ++i) {
		if (chan[i] == "Select") 
{
			channelSelected = Array.concat(channelSelected, false);
		} else
 {
			channelSelected = Array.concat(channelSelected, true);
		}
	}
	return channelSelected;
}
	
/*
 * Remove slice from stack.
 */
function Slice_Remove(chan) {
	for (i=channels; i>0; i--) {
		Stack.setChannel(i);
		if (chan[i-1] == 0)
			run("Delete Slice", "delete=channel");
	}
}

/*
 * Creates a selection from information stored in a file.
 */
function makeSelectionFromFile(file) {
	selectionParameters = File.openAsString(file);
	selectionParameters = split(selection, ",");
	for (i=0; i<selectionParameters.length; ++i) { 
		selectionParameters[i] = parseInt(selectionParameters[i]);
	}
	width = selectionParameters[0];
	height = selectionParameters[1];
	channel = selectionParameters[2];
	x = selectionParameters[3];
	y = selectionParameters[4];
	makeRectangle(x, y, width, height);
	Stack.setChannel(channel);
}

/*
 * Stores rectangular selection in file.
 */
function saveSelectionToFile(file) {
	getSelectionBounds(x, y, width, height);
	Stack.getPosition(channel, VOID, VOID);
	
	// Save to file
	File.saveString("&width,&height,&channel,&x,&y", file);
}

function convertWavelengthToChannelNumber(wavelength, exicationWavelengths) {
	// Convert from wavelength to channel number
	i = 0;
	while (wavelength != excitationWavelengths[i] && i < excitationWavelengths.length) {
		i++;
	}
	return i+1; // index shift
}
