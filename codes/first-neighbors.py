

import numpy as np

# Atom in pixel data structure
class PixelAtom:
  def __init__(self, atom_id=0, coord=None):
    self.atom_id = atom_id                                           # the atom ID
    self.coord = np.zeros(3) if coord is None else np.array(coord)   # the atom coordinates on x, y and z

# Pixel data structure
class Pixel:
  def __init__(self, pid=0, p_xyz=None, tested=False, patoms=0, pix_atoms=None, neighbors=0):
    self.pid = pid                                                  # the pixel number
    self.p_xyz = np.zeros(3) if p_xyz is None else np.array(p_xyz)  # the pixel coordinates in the grid
    self.tested = tested                                            # was the pixel checked already
    self.patoms = patoms                                            # number of atom(s) in pixel
    self.pix_atoms = [] if pix_atoms is None else pix_atoms         # list of atom(s) in the pixel
    self.neighbors = neighbors                                      # number of neighbors for pixel
    self.pixel_neighbors = np.zeros(27, dtype=int)                  # the list of neighbor pixels, maximum 27

# Pixel grid data structure
class PixelGrid:
  def __init__(self, pixels=0, n_pix=None, n_xy=0, pixel_list=None):
    self.pixels = pixels                                            # total number of pixels in the grid
    self.n_pix = np.zeros(3, dtype=int) if n_pix is None else np.array(n_pix) # pixel(s) on each axis
    self.n_xy = n_xy                                                # number of pixels in the plan xy
    self.pixel_list = [] if pixel_list is None else pixel_list      # pointer to the pixels, to be allocated

# Bond distance data structure
class Distance:
  def __init__(self, length=0.0, Rij=None):
    self.length = length                                            # the distance in |\AA| squared
    self.Rij = np.zeros(3) if Rij is None else np.array(Rij)        # vector components of x, y and z

# Model description
atoms = 0                         # the total number of atom(s)
c_coord = None                    # list of Cartesian coordinates: c_coord[atoms][3]
cutoff = 0.0                      # the cutoff to define atomic bond(s)
cutoff_squared = 0.0              # squared value for the cutoff

# Model box description
l_params = np.zeros(3)            # lattice a, b and c
cart_to_frac = np.zeros((3, 3))   # Cartesian to fractional coordinates matrix
frac_to_cart = np.zeros((3, 3))   # fractional to Cartesian coordinates matrix


# Adjust, if needed, shift to search for pixel neighbor(s) using PBC
# - grid pixel_grid          : the pixel grid
# - int pixel_coord[3]       : the pixel coordinates in the grid
# - int pbc_shift[3][3][3]   : the shift, correction, to be calculated
def set_pbc_shift(pixel_grid : PixelGrid, pixel_coord : np.ndarray, pbc_shift : np.ndarray):
  # Initialize pbc_shift to zero
  for x_pos in range(3):
    for y_pos in range(3):
      for z_pos in range(3):
        pbc_shift[x_pos][y_pos][z_pos] = 0        # initialization without any shift

  if pixel_coord[0] == 0:                         # pixel position on 'x' is min
    for y_pos in range(3):
      for z_pos in range(3):
        pbc_shift[0][y_pos][z_pos] = pixel_grid.n_pix[0]

  elif pixel_coord[0] == pixel_grid.n_pix[0] - 1: # pixel position on 'x' is max
    for y_pos in range(3):
      for z_pos in range(3):
        pbc_shift[2][y_pos][z_pos] = -pixel_grid.n_pix[0]

  if pixel_coord[1] == 0:                         # pixel position on 'y' is min
    for x_pos in range(3):
      for z_pos in range(3):
        pbc_shift[x_pos][0][z_pos] += pixel_grid.n_xy

  elif pixel_coord[1] == pixel_grid.n_pix[1] - 1: # pixel position on 'y' is max
    for x_pos in range(3):
      for z_pos in range(3):
        pbc_shift[x_pos][2][z_pos] -= pixel_grid.n_xy

  if pixel_coord[2] == 0:                         # pixel position on 'z' is min
    for x_pos in range(3):
      for y_pos in range(3):
        pbc_shift[x_pos][y_pos][0] += pixel_grid.pixels

  elif pixel_coord[2] == pixel_grid.n_pix[2] - 1: # pixel position on 'z' is max
    for x_pos in range(3):
      for y_pos in range(3):
        pbc_shift[x_pos][y_pos][2] -= pixel_grid.pixels


def add_atom_to_pixel (the_pixel : Pixel, pixel_coord : , atom_id : int, atom_coord)

  if ! the_pixel.patoms:
    for axis in range(3):
      the_pixel.p_xyz[axis] = pixel_coord[axis]
  
  else

  the_pixel.pix_atoms[patoms].atom_id = atom_id
  for axis in range(3):
    the_pixel.pix_atoms[patoms].coord[axis] = atom_coord[axis]

  the_pixel.patoms += 1
 



# Preparation of the pixel grid
# - bool use_pbc : flag to set if PBC are used or not
def prepare_pixel_grid(use_pbc : bool):
  grid = PixelGrid()                                 # Create a new pixel grid
  cmin = [float(|\textquotesingle|inf|\textquotesingle|)] * 3                          # Initialize to infinity
  cmax = [-float(|\textquotesingle|inf|\textquotesingle|)] * 3                         # Initialize to negative infinity
  pixel_pos = np.zeros(3, dtype=int)
  
  if not use_pbc:                                    # Without periodic boundary conditions
    for axis in range(3):
      cmin[axis] = cmax[axis] = c_coord[0][axis]
    for aid in range(1, atoms):                      # For all atoms
      for axis in range(3):                          # For x, y and z
        cmin[axis] = min(cmin[axis], c_coord[aid][axis])
        cmax[axis] = max(cmax[axis], c_coord[aid][axis])
    for axis in range(3):                            # For x, y and z
      grid.n_pix[axis] = int((cmax[axis] - cmin[axis]) / cutoff) + 1  # Number of pixels on axis 'axis'
  else:                                              # Using periodic boundary conditions
    for axis in range(3):                            # For x, y and z
      grid.n_pix[axis] = int(l_params[axis] / cutoff) + 1  # Number of pixels on axis 'axis'
  
  for axis in range(3):                              # For x, y and z
    # Correction if the number of pixels on 'axis' is too small
    grid.n_pix[axis] = 1 if grid.n_pix][axis] < 4 else grid.n_pix[axis]
  
  grid.n_xy = grid.p_pix[0] * grid.n_pix[1] # Number of pixels on the plan 'xy'
  grid.pixels = grid.n_xy * grid.p_pix[2]  # Total number of pixels in the grid
  
  # User defined function to allocate the memory to store the pixel information for the grid
  grid.pixel_list = allocate_pixel_data(grid.pixels)
  
  if not use_pbc:                                    # Without periodic boundary conditions
    for aid in range(atoms):                         # For all atoms
      for axis in range(3):                          # For x, y and z
        pixel_pos[axis] = int((c_coord[aid][axis] - cmin[axis]) / cutoff)
      pixel_num = pixel_pos[0] + pixel_pos[1] * grid.n_pix[0] + pixel_pos[2] * grid.n_xy + 1
      # User defined function to:
      # - Add atom 'aid' with coordinates 'c_coord[aid]' to pixel 'pixel_number'
      # - Increment the number of atom(s) in pixel 'pixel_number'
      # - If needed (for the first atom) set pixel coordinates in the grid to 'pixel_pos'
      add_atom_to_pixel(grid, pixel_num, pixel_pos, aid, c_coord[aid])
  else:  # Using periodic boundary conditions
    for aid in range(atoms):                         # For all atoms
      # with 'matrix_multiplication' a user defined function to perform the operation
      f_coord = matrix_multiplication(cart_to_frac, c_coord[aid])
      for axis in range(3):                          # For x, y and z
        f_coord[axis] = f_coord[axis] - np.floor(f_coord[axis])
        pixel_pos[axis] = int(f_coord[axis] * grid.n_pix[axis])
      pixel_num = pixel_pos[0] + pixel_pos[1] * grid.n_pix[0] + pixel_pos[2] * grid.n_xy + 1
      add_atom_to_pixel(grid, pixel_num, pixel_pos, aid, f_coord)  # User defined function (see above)
  
  return grid



# Finding neighbor pixels for pixel in the grid
# - bool use_pbc : flag to set if PBC are used or not
# - grid * the_grid : pointer to the pixel grid
# - pixel * the_pix : pointer to the pixel with neighbors to be found
def find_pixel_neighbors(use_pbc : bool, the_grid : PixelGrid, the_pix : Pixel):
  boundary = False                            # is pixel on the boundary of the grid
  keep_neighbor = True                        # keep or not neighbor during analysis
  l_start = [0, 0, 0]                         # loop iterators starting value
  l_end = [3, 3, 3]                           # loop iterators ending value
  pmod = [-1, 0, 1]                           # position modifiers
  pbc_shift = np.zeros((3, 3, 3), dtype=int)  # shift for pixel neighbor number due to PBC

  # Check if PBC are used
  if use_pbc:
    set_pbc_shift(the_grid, the_pix.p_xyz, pbc_shift)
  else:
    for axis in range(3):
      if the_pix.p_xyz[axis] == 0 or the_pix.p_xyz[axis] == the_grid.n_pix[axis] - 1:
        boundary = True

  # Adjust the loop start and end based on the grid dimensions
  for axis in range(3):
    if the_grid.n_pix[axis] == 1:
      l_start[axis] = 1
      l_end[axis] = 2

  nnp = 0  # number of neighbors
  for xpos in range(l_start[0], l_end[0]):
    for ypos in range(l_start[1], l_end[1]):
      for zpos in range(l_start[2], l_end[2]):
        keep_neighbor = True

        if not use_pbc and boundary:
          if the_pix.p_xyz[0] == 0 and xpos == 0:
            keep_neighbor = False
          elif the_pix.p_xyz[0] == the_grid.n_pix[0] and xpos == 2:
            keep_neighbor = False
          elif the_pix.p_xyz[1] == 0 and ypos == 0:
            keep_neighbor = False
          elif the_pix.p_xyz[1] == the_grid.n_pix[1] and ypos == 2:
            keep_neighbor = False
          elif the_pix.p_xyz[2] == 0 and zpos == 0:
            keep_neighbor = False
          elif the_pix.p_xyz[2] == the_grid.n_pix[2] and zpos == 2:
            keep_neighbor = False

        if keep_neighbor:
          # Calculate the neighbor id
          nid = the_pix.pid + pmod[xpos] + pmod[ypos] * the_grid.n_pix[0] + pmod[zpos] * the_grid.n_xy
          if use_pbc:
            nid += pbc_shift[xpos][ypos][zpos]
          the_pix.neighbor_list[nnp] = nid
          nnp += 1

  the_pix.neighbors = nnp



# Evaluating the interatomic distance between 2 pixel atoms
# - bool use_pbc : flag to set if PBC are used or not
# - pixel_atom * at_i : pointer to first pixel atom
# - pixel_atom * at_j : pointer to second pixel atom
def evaluate_distance(use_pbc : bool, at_i : PixelAtom, at_j : PixelAtom):
  dist = Distance()         # Placeholder for the distance data structure
  Rij = np.zeros(3)  # Initialize the distance vector
  # Calculating the distance components between atoms
  for axis in range(3):
    Rij[axis] = at_i.coord[axis] - at_j.coord[axis]

  if use_pbc:
    # Pixel atom's coordinates are in corrected fractional format
    for axis in range(3):
      # Absolute value in float format
      u = abs(Rij[axis])
      v = min(u, 1.0 - u)
      # Proper value, with proper sign
      Rij[axis] = (Rij[axis] / u) * v
    
    # Transform back to Cartesian coordinates
    # matrix_multiplication is assumed to be defined elsewhere
    Rij = matrix_multiplication(frac_to_cart, Rij)

  # Calculating the squared distance (no square root for efficiency)
  dist.length = np.sum(Rij ** 2)
  dist.Rij = Rij  # Store the distance vector components

  # Returning the 'distance' data structure that contains:
  # - the squared value for Dij
  # - the components of the distance vector on x, y, and z
  return dist



def pixel_search_for_neighbors(use_pbc : bool):
  # Pointer to the pixel grid for analysis
  all_pixels = prepare_pixel_grid(use_pbc)
  
  # For all pixels in the grid
  for pix in range(all_pixels.pixels):
    # Setting pix_i as pointer to pixel number pix
    pix_i = all_pixels.pixel_list[pix]
    
    # If pixel pix_i contains atom(s)
    if pix_i.patoms:
      # Search for neighbor pixels
      find_pixel_neighbors(use_pbc, all_pixels, pix_i)
      
      # Testing all pix_i neighbor pixels
      for pid in range(pix_i.neighbors):
        pjx = pix_i.pixel_neighbors[pid]
        # Setting pix_j as pointer to pixel number pjx
        pix_j = all_pixels.pixel_list[pjx]
        
        # Checking pixel pix_j if it:
        # - contains atom(s)
        # - was not tested, otherwise the analysis would have been performed already
        if pix_j.patoms and not pix_j.tested:
          # If pix_i and pix_j are the same, only test pair of different atoms
          end = 0 if pjx != pix else 1
          
          # For all atom(s) in pix_i
          for aid in range(pix_i.patoms - end):
            # Set pointer to the first atom to test
            at_i = pix_i.pix_atoms[aid]
            start = 0 if pjx != pix else aid + 1
            
            # For all atom(s) in pix_j
            for bid in range(start, pix_j.patoms):
              # Set pointer to the second atom to test
              at_j = pix_j.pix_atoms[bid]
              
              # Evaluate interatomic distance
              Dij = evaluate_distance(use_pbc, at_i, at_j)
              
              if Dij.length < cutoff_squared:
                # This is a bond!
                pass
                
      # Store that pixel pix_i was tested
      pix_i.tested = True
