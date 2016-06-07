//The lines for formatting the zoom and speed are 206-214

close("*");
print("\\Clear");
run("Set Measurements...", "area mean modal min center bounding median redirect=None decimal=9");

setBatchMode(true);

//Ask user to choose the input and output directories
directory = getDirectory("Choose input directory");
rawFileList = getFileList(directory);
outputDirectory = getDirectory("Choose output directory");
run("Bio-Formats Macro Extensions");

manualOverride = getBoolean("Would you like to manually check each bleach selection area?");

//Count how many files in the directory end in *.lif
lifCount = 0;
for(a=0; a<rawFileList.length; a++){
	if(endsWith(rawFileList[a], ".lif")){
		lifCount = lifCount + 1;
	}
}

//Build a new array of just the *.lif files in the directory
fileList = newArray(lifCount);
lifCount = 0;
for(a=0; a<rawFileList.length; a++){
	if(endsWith(rawFileList[a], ".lif")){
		fileList[lifCount] = rawFileList[a];
		lifCount = lifCount + 1;
	}
}
print("Processing the following files:");
Array.print(fileList);

//Count the number of images there are to be processed within the directory
nSamples = 0;
for (i=0; i<fileList.length; i++) {
	file = directory + fileList[i];
	Ext.setId(file);
	
	//Measure number of series in file
	Ext.getSeriesCount(nStacks);
	nSamples = nSamples + nStacks;
}

//Create an array for saving all sample sames within the directory
sampleNames = newArray(nSamples);
	
//Build an array of all sample names within the directory
sampleCounter = 0;
excludeCounter = 0;
for (i=0; i<fileList.length; i++) {
	file = directory + fileList[i];
	Ext.setId(file);
	
	//Measure number of series in file
	Ext.getSeriesCount(nStacks);

	//Open all stacks from set of lif files, one stack at a time
	for(j=0; j<nStacks; j++) {	

		//Only open the image if it has one / (i.e. part of collection, but not subcollection), and ends with a series number (not a, b, merged, etc.)
		Ext.setSeries(j);
		Ext.getSeriesName(seriesTitle);

		//Make sure the series is the pre-bleach part
		if(matches(seriesTitle, ".*/.*/FRAP Pre Series.*$")){
			seriesTitle = replace(seriesTitle, fileList[i] + " - ", "");
			sampleNames[sampleCounter] = replace(seriesTitle, "/.*/.*", "");
		}
		//If this series name is the post bleach portion, then remove name (identical experiment)
		else if (matches(seriesTitle, ".*/.*/FRAP Pb1 Series.*$")){
			sampleNames[sampleCounter] = "N/A";
		}	
		//If this series name is not formatted correctly, exclude it from the list
		else{
			sampleNames[sampleCounter] = "N/A";
			excludeCounter = excludeCounter + 1;
		}
			
		//Count up one sample index
		sampleCounter = sampleCounter + 1;
	}
}

print(excludeCounter + " images/stacks were excluded from analysis due to incorrect naming format.");

//Record the number of series within the same sample:
sampleNamesCopy = Array.copy(sampleNames);
seriesNumbers = newArray(sampleNames.length);
Array.fill(seriesNumbers, 0);
//Search for a unique sample name
for(a=0; a<sampleNames.length; a++){
	seriesIDcounter = 1;
	
	//If found, replace name with "N/A" and record series number as 1
	if(sampleNamesCopy[a] != "N/A"){
		currentName = sampleNamesCopy[a];
		sampleNamesCopy[a] = "N/A";
		seriesNumbers[a] = seriesIDcounter;
		seriesIDcounter = seriesIDcounter + 1;

		//Search for all identical sample names and replace with N/A and record series number accordingly
		for(b=a; b<sampleNames.length; b++){
			if(sampleNamesCopy[b] == currentName){
				sampleNamesCopy[b] = "N/A";
				seriesNumbers[b] = seriesIDcounter;
				seriesIDcounter = seriesIDcounter + 1;
			}
		}
	}
}

//Open all data contained within the directory
sampleCounter = 0;

print("");
print("Series_Name,Sample,Bleach_Correlation,Surround_Correlation,delta_Int,Area,Bleach_X,Bleach_Y,Bleach_Width,Bleach_Height,Area_Set_By");

for (i=0; i<fileList.length; i++) {
	file = directory + fileList[i];
	Ext.setId(file);

	//Measure number of series in file
	Ext.getSeriesCount(nStacks);
	
	//Open corresponding pair of lif files
	for(j=0; j<nStacks; j++) {	

		//Only open the image if it was formatted correctly (i.e. not listed as "N/A")
		if(sampleNames[sampleCounter] != "N/A"){

			seriesName = sampleNames[sampleCounter] + " - Series " + seriesNumbers[sampleCounter];
			showProgress(sampleCounter/nSamples);

			//Open pre-bleach and post-bleach pair of stacks
			run("Bio-Formats Importer", "open=file color_mode=Default view=Hyperstack stack_order=XYCZT series_"+d2s(j+1,0));
			preBleach = getTitle();
			run("Bio-Formats Importer", "open=file color_mode=Default view=Hyperstack stack_order=XYCZT series_"+d2s(j+2,0));
			postBleach = getTitle();

			//Set properties to pixel units
			run("Properties...", "unit=pixel pixel_width=1 pixel_height=1 voxel_depth=1");
	
			//Combine the pre-bleach stack with the post bleach stack
			run("Concatenate...", "  title=[Concatenated Stacks.tif] image1=[" + preBleach + "] image2=[" + postBleach + "]");
			saveAs("Tiff", outputDirectory + "Concatenated Stacks.tif");
			
			//Split the channels to allow binning og time points (hyperstacks can bin time points)
			run("Split Channels");
			selectWindow("C2-Concatenated Stacks.tif");
			run("Bin...", "x=1 y=1 z=3 bin=Average");
			
			//Create a stdev z-projection to find the bleach region
			run("Gaussian Blur...", "sigma=10 stack");
			run("Z Project...", "start=2 projection=[Average Intensity]");
			selectWindow("C2-Concatenated Stacks.tif");	
			run("Stack to Images");
			imageCalculator("Subtract create", "C2-Concatenated-0001","AVG_C2-Concatenated Stacks.tif");
			setAutoThreshold("Otsu dark");
			run("Analyze Particles...", "size=0-Infinity pixel display clear");

			//look for the area the is most likely the bleach area
			areaScore = 0;
			bleachSpot = -1;
			for(a=0; a<nResults; a++){
				bleachArea = getResult("Area",a);
				bleachMean = getResult("Mean",a);
				if(bleachMean * bleachArea > areaScore){
					areaScore = bleachMean * bleachArea;
					bleachSpot = a;
				}
			}
			
			//Close all images and re-open them
			close("*");

			//If the manual override is on, turn off batch mode to allow user to see sample
			if(manualOverride){
				setBatchMode(false);
			}
			
			open(outputDirectory + "Concatenated Stacks.tif");
			
			//If a candidate bleach region was found, analyze it
			if(bleachSpot != -1){
				//Create a selection based off of the bleach area
				makeRectangle(getResult("BX",bleachSpot), getResult("BY",bleachSpot), getResult("Width",bleachSpot), getResult("Height",bleachSpot));

				//Create arrays for measuring the correlation
				CFPinBleach = newArray(15);
				CFPoutBleach = newArray(15);
				YFPinBleach = newArray(15);
				YFPoutBleach = newArray(15);
				
				//Count the number of frames in the stack
				Stack.getDimensions(dummy, dummy, dummy, dummy, nFrames) 

				//If the manual override option was selected, then allow user to view and then override selection area
				if(manualOverride){
					//Allow user to see the YFP selection area and verify the seletion
					Stack.setChannel(2);
					run("Enhance Contrast", "saturated=0.1");
					run("In [+]");
					run("In [+]");
					run("In [+]");
					wait(500);
					
					for(frame = 1; frame<=nFrames; frame++){
						Stack.setFrame(frame);
						wait(100);
					}
	
					//Ask the user to verify that the selection is correct
					waitForUser("Adjust the selection area if necessary");
	
					//Get the selection bounds
					getSelectionBounds(x, y, width, height);

					//If the selection was changed, then record the new selection area
					if(x != getResult("BX",bleachSpot) || y != getResult("BY",bleachSpot) || width != getResult("Width",bleachSpot) || height != getResult("Height",bleachSpot)){
						setResult("BX",bleachSpot,x);
						setResult("BY",bleachSpot,y);
						setResult("Width",bleachSpot,width);
						setResult("Height",bleachSpot,height);
	
						//record that the user set the bleach area
						bleachSetBy = "User";
					}
					//Otherwise, keep the computer area, and record that the computer set the area
					else{
						bleachSetBy = "Computer";
					}

					//Turn batch mode back on for the remainder of the analysis
					setBatchMode(true);
				}
				//If no override, record that data was acquired by computer
				else{
					bleachSetBy = "Computer";
				}

				//If there is a selection, continue with the analysis
				if(selectionType() == 0){
					//Collect all intensity measurements
					for(selection = 0; selection <=1; selection++){
						//Invert the selection
						run("Make Inverse");
						for(channel=0; channel<=1; channel++){
							//Set the channel accordingly
							Stack.setChannel(channel+1);
							for(frame=0; frame<nFrames; frame++){
								Stack.setFrame(frame+1)
								getStatistics(dummy, mean);
								
								//Store the mean intensity measurement in the correct array
								if(selection && channel) CFPinBleach[frame] = mean;
								if(selection && !channel) YFPinBleach[frame] = mean;
								if(!selection && channel) CFPoutBleach[frame] = mean;
								if(!selection && !channel) YFPoutBleach[frame] = mean;
							}
						}
					}
					
					//Get the means of each array
					Array.getStatistics(CFPinBleach, dummy, dummy, CFPinMean, dummy);
					Array.getStatistics(YFPinBleach, dummy, dummy, YFPinMean, dummy);
					Array.getStatistics(CFPoutBleach, dummy, dummy, CFPoutMean, dummy);
					Array.getStatistics(YFPoutBleach, dummy, dummy, YFPoutMean, dummy);
					
					//Measure the bleach area Pearson product-moment correlation coefficient
					corrNum = 0;
					corrDenom1 = 0;
					corrDenom2 = 0;
					
					for(a=0; a<nFrames; a++){
						corrNum = corrNum + (CFPinBleach[a] - CFPinMean)*(YFPinBleach[a] - YFPinMean);
						corrDenom1 = corrDenom1 + (CFPinBleach[a] - CFPinMean)*(CFPinBleach[a] - CFPinMean);
						corrDenom2 = corrDenom2 + (YFPinBleach[a] - YFPinMean)*(YFPinBleach[a] - YFPinMean);
					}
					
					//Calculate the correlation coefficient
					bleachCorr = corrNum / (sqrt(corrDenom1)*sqrt(corrDenom2));
					
					//Measure the non-bleach area Pearson product-moment correlation coefficient
					corrNum = 0;
					corrDenom1 = 0;
					corrDenom2 = 0;
					
					for(a=0; a<nFrames; a++){
						corrNum = corrNum + (CFPoutBleach[a] - CFPoutMean)*(YFPoutBleach[a] - YFPoutMean);
						corrDenom1 = corrDenom1 + (CFPoutBleach[a] - CFPoutMean)*(CFPoutBleach[a] - CFPoutMean);
						corrDenom2 = corrDenom2 + (YFPoutBleach[a] - YFPoutMean)*(YFPoutBleach[a] - YFPoutMean);
					}
					
					//Calculate the correlation coefficient
					outerCorr = corrNum / (sqrt(corrDenom1)*sqrt(corrDenom2));
	
					//Print results
					print(preBleach + "," + seriesName + "," + bleachCorr + "," + outerCorr + "," + getResult("Mean",bleachSpot) + "," + getResult("Area",bleachSpot) + "," + getResult("BX",bleachSpot) + "," + getResult("BY",bleachSpot)+ "," + getResult("Width",bleachSpot) + "," + getResult("Height",bleachSpot) + "," + bleachSetBy);
	
				}
				//If there was no selection then record the data accordingly
				else{
					print(preBleach + "," + seriesName + ",N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,User");
				}
			}
			//If no bleach spot was found, report N/A
			else{
				print(preBleach + "," + seriesName + ",N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,Computer");
			}
	
				//Increment the sample counter by 2
				sampleCounter = sampleCounter + 2;
	
				//Increment the loop counter one extra step
				j = j + 1;
			
				close("*");
		}
		else{
			//Increment the sample counter by 1
			sampleCounter = sampleCounter + 1;			
		}
	}
}

//Delete the stack file and save the results
dummy = File.delete(outputDirectory + "Concatenated Stacks.tif");

selectWindow("Log");
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
logName = "" + year + "-" + month + "-" + dayOfMonth + " FRET analysis at " + hour + " hours " + minute + " minutes.csv";
saveAs("Text", outputDirectory + logName);
run("Close");

selectWindow("Results");
run("Close");

setBatchMode(false);