# Import modules
import os
import sys
import arcgisscripting
import arcpy
import traceback

# Create the geoprocessor object.
gp = arcgisscripting.create()

# Set the necessary product code.
gp.SetProduct("ArcInfo")

# Check out necessary licenses.
gp.CheckOutExtension("Spatial")
gp.CheckOutExtension("Network")

# Overwrite output.
gp.overwriteoutput = True


#############################################


def fold_creater(workspace, data_list):
    ''' Creates folders for each output file type'''
    os.chdir(workspace)
    # Create folders for each of the shapefile types
    for item in data_list:
        # Check if folder exists
        if os.path.isdir(item):
            arcpy.AddMessage("The Folder " + item + " already exists!")

        # Create new folder if it doesn't exist
        else:
            arcpy.AddMessage("The folder " + item + " has been created")
            arcpy.CreateFolder_management(workspace, item)


def merge_shp(data_list):
    ''' Merges all of the shapefiles for every initial folders'''
    paths = []
    for data in data_list:
        # Get the list of files
        current_ws = os.path.join(orig_ws + "\\" + data)
        paths.append(current_ws)

    for folder in paths:
        to_merge = []
        # Extract the name of the folder
        rawname = folder.split("\\")[-1]

        # Extract all shapefiles
        for sfile in os.listdir(folder):
            if sfile.endswith(".shp"):
                to_merge.append(sfile)

        # Change the workspace to the new folder
        basepath = (os.path.join(orig_ws, rawname))
        os.chdir(basepath)

        # allow for unique output name of the file
        outname = os.path.join(output_ws, rawname + "\\" + rawname + '_m' + '.shp')
        unique = arcpy.CreateUniqueName(outname, basepath)
        arcpy.Merge_management(to_merge, unique)


def merge_tbl(table_list):
    ''' Merges all of the tables for every initial geodatabases'''
    for gdb in table_list:
        gdb_ws = (os.path.join(orig_ws, gdb))
        arcpy.env.workspace = gdb_ws

        # List all tables in each of the original gdbs
        tables = arcpy.ListTables()

        # Get the name of the folder
        gdbname = gdb.split(".")[0]

        # Create the path for the output folder
        outpath = os.path.join(output_ws, gdbname)

        # Create a unique name for the output file
        unique = arcpy.CreateUniqueName(gdbname + '_m' +'.dbf', outpath)

        # Merge all input files and save to the output folder
        arcpy.Merge_management(tables, unique)
#############################################


def main(output_ws, orig_ws):
    # List of the files created with the VMT
    data_list = ['elevationpoints', 'longest_riverlines', 'midcrosslines', 'rivercrosslines', 'valleyfloors',
                 'watershedmidlines']
    table_list = ['elevtables.gdb', 'midcrosstbls.gdb', 'vfwtables.gdb']

    # Ensure the output folder has all needed shapefile folders
    fold_creater(output_ws, data_list)
    arcpy.AddMessage("The shapefile folders have been created ...")

    # Ensure the output folder has all needed table folders

    # Remove all extenstions into new list
    cleantbl_list = []
    for item in table_list:
        cleantbl_list.append(item.split('.')[0])

    # Create folders without extensions
    fold_creater(output_ws, cleantbl_list)
    arcpy.AddMessage("The table folders have been created ...")

    # Merge all shapefiles
    merge_shp(data_list)
    arcpy.AddMessage("The shapefiles have been merged ...")

    # Merge all tables
    merge_tbl(table_list)
    arcpy.AddMessage("The tables have been merged ...")


    os.chdir(orig_ws)
#############################################


if __name__ == '__main__':
    try:
        # Get the input parameters
        orig_ws = arcpy.GetParameterAsText(0)  # Workspace for the VMT data
        output_ws = arcpy.GetParameterAsText(1)  # Output folder for merged files
        gp.scratchworkspace = arcpy.GetParameterAsText(2)  # Folder for temp files

        # Set the workspace to the VMT data
        arcpy.env.workspace = orig_ws

        # Generate merged files
        main(output_ws, orig_ws)

    except:
        gp.addmessage(traceback.format_exc())
        del gp
        sys.exit()
