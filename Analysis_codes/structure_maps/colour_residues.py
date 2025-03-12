"""
colour_residues.py

Pymol script to map normalised relative errors to the 3D protein structure in the form of a heat map.
Uses a dictionary of PDB IDs mapped to a list of (residue number, relative error) tuples saved as a .pkl file (created using map_devs.py) for each NMR observable.

run colour_residues.py

"""



from pymol import cmd
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import sys
import pickle
# run ../colour_residues.py

prot = '1ozi' # Change 
norm_lim = 0.5 # Change to 1 if necessary

def color_residues_by_error1(residue_errors, prot):
    """
    Colors residues based on error values using a color map.
    
    Parameters:
    residue_errors (list of tuples): List of (residue_number, error_value) tuples.
    """
    # Normalize error values to range [0, 1]
    errors = [abs(error) for _, error in residue_errors]
    min_error = min(errors)
    max_error = max(errors)
    normalized_errors = [(res - min_error) / (max_error - min_error) for res in errors]
    
    # Create a color map
    cmd.color("blue", prot)
    cmap = plt.get_cmap('coolwarm')
    
    # Apply colors to residues
    for (residue, error), norm_error in zip(residue_errors, normalized_errors):
        color = cmap(norm_error)
        color_name = f"color_{residue}"
        cmd.set_color(color_name, [color[0], color[1], color[2]])
        cmd.color(color_name, f"resi {residue}")
    
    fig, ax = plt.subplots(figsize=(6, 1))
    fig.subplots_adjust(bottom=0.5)
    
    norm = plt.Normalize(min_error, max_error)
    cb1 = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=ax, orientation='horizontal')
    cb1.set_label('Normalised Relative Errors')
    
    plt.savefig('color_bar.png')
    plt.close()

def color_residues_by_error(residue_errors, vmin=None, vmax=None, prot):
    """
    Colors residues based on error values using a color map.
    
    Parameters:
    residue_errors (list of tuples): List of (residue_number, error_value) tuples.
    vmin (float, optional): Minimum value for the color scale. Defaults to None.
    vmax (float, optional): Maximum value for the color scale. Defaults to None.
    """
    # Normalize error values to range [0, 1]
    errors = [abs(error) for _, error in residue_errors]
    min_error = min(errors) if vmin is None else vmin
    max_error = max(errors) if vmax is None else vmax
    normalized_errors = [(res - min_error) / (max_error - min_error) for res in errors]
    
    # Create a color map
    cmd.color("blue", prot)
    cmap = plt.get_cmap('RdYlGn_r')
    
    # Apply colors to residues
    for (residue, error), norm_error in zip(residue_errors, normalized_errors):
        color = cmap(norm_error)
        color_name = f"color_{residue}"
        cmd.set_color(color_name, [color[0], color[1], color[2]])
        cmd.color(color_name, f"resi {residue}")
    
    fig, ax = plt.subplots(figsize=(6, 1))
    fig.subplots_adjust(bottom=0.5)
    
    norm = plt.Normalize(vmin=min_error, vmax=max_error)
    cb1 = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=ax, orientation='horizontal')
    cb1.set_label('Absolute values of Relative Errors')
    
    # plt.savefig('../rel_err_shifts/color_bar.png')

def save_image_with_color_bar(pymol_image, color_bar_image, output_image):
    """
    Combine the PyMOL image and the color bar image into one image and save it.
    
    Parameters:
    pymol_image (str): The filename of the PyMOL image.
    color_bar_image (str): The filename of the color bar image.
    output_image (str): The filename of the combined output image.
    """
    pymol_img = Image.open(pymol_image)
    color_bar_img = Image.open(color_bar_image)
    
    # Create a new image with enough space for both images
    total_width = pymol_img.width
    total_height = pymol_img.height + color_bar_img.height
    
    combined_img = Image.new('RGB', (total_width, total_height))
    
    # Paste the images into the combined image
    combined_img.paste(pymol_img, (0, 0))
    combined_img.paste(color_bar_img, (0, pymol_img.height))
    
    # Save the combined image
    combined_img.save(output_image)

def save_image(filename):
    """
    Save the current PyMOL session as an image.
    
    Parameters:
    filename (str): The name of the file to save the image as.
    """
    cmd.png(filename)

# Main code to load relative errors and colour residues in pymol

with open('../../rel_err_shifts_H.pkl', 'rb') as f:
    residue_errors = pickle.load(f)

cmd.bg_color("white")
# color_residues_by_error(residue_errors[prot], 0, 1)
# save_image("../rel_err_shifts/"+prot+"_colored_residues_norm1.png")
color_residues_by_error(residue_errors[prot], 0, norm_lim, prot)
# save_image("../structures_charmm/"+prot+"_colored_residues.png")
# save_image_with_color_bar("../structures_charmm/"+prot+"_colored_residues.png", "../structures_charmm/color_bar.png", "../structures_charmm/"+prot+"_combined_image.png")


