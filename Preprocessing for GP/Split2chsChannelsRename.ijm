// Batch process .nd2 files in a selected directory
dir = getDirectory("Choose a Directory Containing .nd2 Files");
list = getFileList(dir);

for (i = 0; i < list.length; i++) {
    if (endsWith(list[i], ".nd2")) {
        open(dir + list[i]); // Open the .nd2 file
        
        // Get filename without extension
        filename = File.nameWithoutExtension;

        // Split channels
        run("Split Channels");

        // Dynamically fetch the names of open windows
        channel1 = "C1-" + list[i];
        channel2 = "C2-" + list[i];
        channel3 = "C3-" + list[i];

        // Rename and save the first channel
        selectWindow(channel1);
        rename(filename + "_ch00");
        saveAs("Tiff", dir + filename + "_ch00.tif");
        close();

        // Rename and save the second channel
        selectWindow(channel2);
        rename(filename + "_ch01");
        saveAs("Tiff", dir + filename + "_ch01.tif");
        close();

        // Close the third channel without saving
        selectWindow(channel3);
        close();
    }
}

print("Batch processing complete.");
