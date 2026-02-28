


#define TRUE 1
#define FALSE 0

// atom in pixel data structure
typedef struct pixel_atom pixel_atom;
struct pixel_atom
{
  int atom_id;             // the atom ID
  float coord[3];          // the atom coordinates on x, y and z
};

// pixel data structure
typedef struct pixel pixel;
struct pixel
{
  int pid;                 // the pixel number
  int p_xyz[3];            // the pixel coordinates in the grid
  bool tested;             // was the pixel checked already
  int patoms;              // number of atom(s) in pixel
  pixel_atom * pix_atoms;  // list of atom(s) in the pixel, to be allocated
  int neighbors;           // number of neighbors for pixel
  int pixel_neighbors[27]; // the list of neighbor pixels, maximum 27
};

// pixel grid data structure
typedef struct pixel_grid pixel_grid;
struct pixel_grid
{
  int pixels;              // total number of pixels in the grid
  int n_pix[3];            // number of pixel(s) on each axis
  int n_xy;                // number of pixels in the plan xy
  pixel * pixel_list;      // pointer to the pixels, to be allocated
};

// bond distance data structure
typedef struct distance distance;
struct distance
{
  float length;            // the distance in \AA\ squared
  float Rij[3];            // vector components of x, y and z
};

/*
  the following are considered to be provided by the user
*/

// model description
int atoms;                 // the total number of atom(s)
float ** c_coord;          // list of Cartesian coordinates: c_coord[atoms][3]
float cutoff;              // the cutoff to define atomic bond(s)
float cutoff_squared;      // squared value for the cutoff

// model box description
float l_params[3];         // lattice a, b and c
float cart_to_frac[3][3];  // Cartesian to fractional coordinates matrix
float frac_to_cart[3][3];  // fractional to Cartesian coordinates matrix


// Adjust, if needed, shift to search for pixel neighbor(s) using PBC
// - pixel_grid * grid           : the pixel grid
// - int pixel_coord[3]         : the pixel coordinates in the grid
// - int pbc_shift[3][3][3]     : the shift, correction, to be calculated
void set_pbc_shift (pixel_grid * grid, int pixel_coord[3], int pbc_shift[3][3][3])
{
  int x_pos, y_pos, z_pos; // loop iterators
  for ( x_pos = 0 ; x_pos < 3 ; x_pos ++ )
  {
    for ( y_pos = 0 ; y_pos < 3 ; y_pos ++ )
    {
      for ( z_pos = 0 ; z_pos < 3 ; z_pos ++ )
      {
        pbc_shift[x_pos][y_pos][z_pos] = 0;           // initialization without any shift
      }
    }
  }
  if ( pixel_coord[0] == 0 )                            // pixel position on 'x' is min
  {
    for ( y_pos = 0 ; y_pos < 3 ; y_pos ++ )
    {
      for ( z_pos = 0 ; z_pos < 3 ; z_pos ++ )
      {
        pbc_shift[0][y_pos][z_pos] = grid->n_pix[0];
      }
    }
  }
  else if ( pixel_coord[0] == grid->n_pix[0] - 1 )      // pixel position on 'x' is max
  {
    for ( y_pos = 0 ; y_pos < 3 ; y_pos ++ )
    {
      for ( z_pos = 0 ; z_pos < 3 ; z_pos ++ )
      {
        pbc_shift[2][y_pos][z_pos] = - grid->n_pix[0];
      }
    }
  }
  if ( pixel_coord[1] == 0 )                            // pixel position on 'y' is min
  {
    for ( x_pos = 0 ; x_pos < 3 ; x_pos ++ )
    {
      for ( z_pos = 0 ; z_pos < 3 ; z_pos ++ )
      {
        pbc_shift[x_pos][0][z_pos] += grid->n_xy;
      }
    }
  }
  else if ( pixel_coord[1] == grid->n_pix[1] - 1 )      // pixel position on 'y' is max
  {
    for ( x_pos = 0 ; x_pos < 3 ; x_pos ++ )
    {
      for ( z_pos = 0 ; z_pos < 3 ; z_pos ++ )
      {
        pbc_shift[x_pos][2][z_pos] -= grid->n_xy;
      }
    }
  }
  if ( pixel_coord[2] == 0 )                            // pixel position on 'z' is min
  {
    for ( x_pos = 0 ; x_pos < 3 ; x_pos ++ )
    {
      for ( y_pos = 0 ; y_pos < 3 ; y_pos ++ )
      {
        pbc_shift[x_pos][y_pos][0] += grid->pixels;
      }
    }
  }
  else if ( pixel_coord[2] == grid->n_pix[2] - 1 )      // pixel position on 'z' is max
  {
    for ( x_pos = 0 ; x_pos < 3 ; x_pos ++ )
    {
      for ( y_pos = 0 ; y_pos < 3 ; y_pos ++ )
      {
        pbc_shift[x_pos][y_pos][2] -= grid->pixels;
      }
    }
  }
}


void add_atom_to_pixel (pixel * the_pixel, int pixel_coord[3], int atom_id, float atom_coord[3])
{
  int axis;
  if (! the_pixel->patoms) 
  {
    // if the pixel do contains any atom yet, then save its coordinates in the grid
    for ( axis = 0 ; axis < 3 ; axis ++ )
    {
      the_pixel->p_xyz[axis] = pixel_coord[axis];
    }
    // allocate the memory to store the first pixel_atom information
    the_pixel->pix_atoms = malloc(sizeof*the_pixel->pix_atoms);
  }
  else
  {
    // otherwise reallocate memory to store the new pixel_atom information
    the_pixel->pix_atoms = realloc(the_pixel->pix_atoms, (the_pixel->patoms+1)*sizeof*the_pixel->pix_atoms);
  }
  the_pixel->pix_atoms[patoms].atom_id = atom_id;
  for ( axis = 0 ; axis < 3 ; axis ++ )
  {
    the_pixel->pix_atoms[patoms].coord[axis] = atom_coord[axis];
  }
  // increment the number of atom(s) in the pixel
  the_pixel->patoms ++;
}


pixel_grid * prepare_pixel_grid (bool use_pbc)
{
  pixel_grid * grid;        // pointer to the pixel grid to create
  int axis;                 // loop iterator axis id (0=x, 1=y, 2=z)
  int aid;                  // loop iterator atom number (0, atoms-1)
  int pixel_num;            // pixel number in the grid
  int pixel_pos[3];         // pixel coordinates in the grid
  float cmin[3], cmax[3];   // float coordinates min, max values
  float f_coord[3];         // float fractional coordinates

  grid = malloc(sizeof*grid);

  if ( ! use_pbc ) // without periodic boundary conditions
  {
    for ( axis = 0 ; axis < 3 ; axis ++ ) cmin[axis] = cmax[axis] = c_coord[0][axis];
    for ( aid = 1 ; aid < atoms ; aid ++ )    // For all atoms
    {
      for ( axis = 0 ; axis < 3 ; axis ++ )   // For x, y and z
      {
        cmin[axis] = min(cmin[axis], c_coord[aid][axis]);
        cmax[axis] = max(cmax[axis], c_coord[aid][axis]);
      }
    }
    for ( axis = 0 ; axis < 3 ; axis ++ )     // For x, y and z
    {
      grid->n_pix[axis] = (int)((cmax[axis] - cmin[axis]) / cutoff) + 1; // Number of pixel(s) on axis 'axis'
    }
  }
  else  // using periodic boundary conditions
  {
    for ( axis = 0 ; axis < 3 ; axis ++ )     // For x, y and z
    {
      grid->n_pix[axis] = (int)(l_params[axis] / cutoff) + 1; // Number of pixel(s) on axis 'axis'
    }
  }
  for ( axis = 0 ; axis < 3 ; axis ++ )       // For x, y and z
  {
    // correction if the number of pixel(s) on 'axis' is too small
    grid->n_pix[axis] = (grid->n_pix[axis] < 4) ? 1 : grid->n_pix[axis];
  }
  grid->n_xy = grid->n_pix[0] * grid->n_pix[1]; // Number of pixels on the plan 'xy'
  grid->pixels = grid->n_xy * grid->n_pix[2];   // Total number of pixels in the grid
  grid->pixel_list = malloc (grid->pixels*sizeof*grid->pixel_list);

  if ( ! use_pbc ) // without periodic boundary conditions
  {
    for ( aid = 0 ; aid < atoms ; aid ++ )  // for all atoms
    {
      for ( axis = 0 ; axis < 3 ; axis ++ )    // for x, y and z
      {
        pixel_pos[axis] = (int)((c_coord[aid][axis] - cmin[axis])/cutoff);
      }
      pixel_num = pixel_pos[0] + pixel_pos[1] * grid->n_pix[0] + pixel_pos[2] * grid->n_xy + 1;
      add_atom_to_pixel (grid, pixel_num, pixel_pos, aid, c_coord[aid]);
    }
  }
  else // using periodic boundary conditions
  {
    for ( aid = 0 ; aid < atoms ; aid ++ )  // for all atoms
    {
      // with 'matrix_multiplication' a user defined function to perform the operation
      f_coord = matrix_multiplication (cart_to_frac, c_coord[aid]);
      for ( axis = 0 ; axis < 3 ; axis ++ )    // for x, y and z
      {
        f_coord[axis] = f_coord[axis] - floorf(f_coord[axis]);
        pixel_pos[axis] = (int)((f_coord[axis] * n_pix[axis]);
      }
      pixel_num = pixel_pos[0] + pixel_pos[1] * grid->n_pix[0] + pixel_pos[2] * grid->n_xy + 1;
      add_atom_to_pixel (& grid->pixel_list[pixel_num], aid, f_coord);
    }
  }
  return grid;
}


// finding neighbor pixels for pixel in the grid
// - bool use_pbc : flag to set if PBC are used or not
// - pixel_grid * the_grid : pointer to the pixel grid
// - pixel * the_pix : pointer to the pixel with neighbors to be found
void find_pixel_neighbors (bool use_pbc, pixel_grid * the_grid, pixel * the_pix)
{
  int axis;                       // loop iterator axis id (1=x , 2=y , 3= z)
  int xpos, ypos, zpos;           // neighbor position on x, y and z
  int l_start[3] = { 0, 0, 0};    // loop iterators starting value
  int l_end[3] = { 3, 3, 3};      // loop iterators ending value
  int pmod[3] = {-1, 0, 1};       // position modifiers
  int nnp;                        // number of neighbors for pixel
  int nid;                        // neighbor id for pixel
  int pbc_shift[3][3][3];         // shift for pixel neighbor number due to PBC
  bool boundary = FALSE;          // is pixel on the boundary of the grid
  bool keep_neighbor = TRUE;      // keep or not neighbor during analysis

  if ( use_pbc )
  {
    set_pbc_shift (the_grid, the_pix->p_xyz, pbc_shift);
  }
  else
  {
    for ( axis = 0 ; axis < 3 ; axis ++ )
    {
      if ( the_pix->p_xyz[axis] == 0 || the_pix->p_xyz[axis] == the_grid->n_pix[axis] - 1 ) boundary = TRUE;
    }
  }
  for ( axis = 0 ; axis < 3 ; axis ++ )
  {
    if ( the_grid->n_pix[axis] == 1 )
    {
      l_start[axis] = 1;
      l_end[axis] = 2;
    }
  }
  nnp = 0;
  for ( xpos = l_start[0] ; xpos < l_end[0] ; xpos ++ )
  {
    for ( ypos = l_start[1] ; ypos < l_end[1] ; ypos ++ )
    {
      for ( zpos = l_start[2] ; zpos < l_end[2] ; zpos ++ )
      {
        keep_neighbor = TRUE;
        if ( ! use_pbc && boundary )
        {
          if (( the_pix->p_xyz[0] == 0 && xpos == 0 ) || ( the_pix->p_xyz[0] == the_grid->n_pix[0] && xpos == 2 ))
          {
            keep_neighbor = FALSE;
          }
          else if (( the_pix->p_xyz[1] == 0 && ypos == 0 ) || ( the_pix->p_xyz[1] == the_grid->n_pix[1] && ypos == 2 ))
          {
            keep_neighbor = FALSE;
          }
          else if (( the_pix->p_xyz[2] == 0 && zpos == 0 ) || ( the_pix->p_xyz[2] == the_grid->n_pix[2] && zpos == 2 ))
          {
            keep_neighbor = FALSE;
          }
        }
        if ( keep_neighbor )
        {
          nid = the_pix->pid + pmod[xpos] + pmod[ypos] * the_grid->n_pix[0] + pmod[zpos] * the_grid->n_xy;
          if ( use_pbc ) nid += pbc_shift[xpos][ypos][zpos];
          the_pix->pixel_neighbors[nnp] = nid;
          nnp ++ ;
        }
      }
    }
  }
  the_pix->neighbors = nnp;
}


// evaluating the interatomic distance between 2 pixel atoms
// - bool use_pbc : flag to set if PBC are used or not
// - pixel_atom * at_i : pointer to first pixel atom
// - pixel_atom * at_j : pointer to second pixel atom
distance evaluate_distance (bool use_pbc, pixel_atom * at_i, pixel_atom * at_j)
{
  int axis;          // axis loop iterator
  float u, v;        // float parameters
  distance dist;     // distance data to store calculation results
  for ( axis = 0 ; axis < 3 ; axis ++ )
  {
    dist.Rij[axis] = at_i->coord[axis] - at_j->coord[axis];
  }
  if ( use_pbc )
  {
    // then the pixel_atom's coordinates are in corrected fractional format
    for ( axis = 0 ; axis < 3 ; axis ++ )
    {
      // absolute value in float format
      u = fabs (dist.Rij[axis]);
      v = min (u, 1.0 - u);
      // proper value, with proper sign
      dist.Rij[axis] = (dist.Rij[axis] / u) * v;
    }
    // transform back to Cartesian coordinates
    // with 'matrix_multiplication' a user defined function to perform the operation
    dist.Rij = matrix_multiplication (frac_to_cart, dist.Rij);
  }
  dist.length = 0.0;
  for ( axis = 0 ; axis < 3 ; axis ++ )
  {
    dist.length += dist.Rij[axis] * dist.Rij[axis];
  }
  // returning the 'distance' data structure that contains:
  // - the squared value for Dij: no time consuming square root calculation !
  // - the components of the distance vector on x, y and z
  return dist;
}


// searching for first neighbor atoms using the grid pixelation/partitioning method
// - bool use_pbc : flag to set if PBC are used or not
void pixel_search_for_neighbors (bool use_pbc)
{
  pixel_grid * all_pixels;   // pointer to the pixel grid for to analyze
  int pix, pjx;              // integer pixel ID numbers
  int aid, bid;              // integer loop atom numbers
  int pid;
  int l_start, l_end;         // integer loop modifier
  pixel * pix_i, * pix_j;     // pointers on pixel data structure
  pixel_atom * at_i, * at_j;  // pointers on pixel_atom data structure
  distance Dij;               // distance data structure

  all_pixels = prepare_pixel_grid (use_pbc);
  // note that the pixel grid 'all_pixels' must be prepared before the following
  // for all pixels in the grid
  for ( pix = 0 ; pix < all_pixels->pixels ; pix ++ )
  {
    // setting 'pix_i' as pointer to pixel number 'pix'
    pix_i = & all_pixels->pixel_list[pix];
    // if pixel 'pix_i' contains atom(s)
    if ( pix_i->patoms )
    {
      // search for neighbbor pixels of 'pix_i'
      find_pixel_neighbors ( use_pbc, all_pixels, pix_i );
      // testing all 'pix_i' neighbor pixels
      for ( pid = 0 ; pid < pix_i->neighbors ; pid ++ )
      {
        pjx = pix_i->pixel_neighbors[pid];
        // setting 'pix_j' as pointer to pixel number 'pjx'
        pix_j = & all_grid->pixel_list[pjx];
        // checking pixel 'pix_j' if it:
        // - contains atom(s)
        // - was not tested, otherwise the analysis would have been performed already
        if ( pix_j->patoms && ! pix_j->tested )
        {
          // if 'pix_i' and 'pix_j' are the same, only test pair of different atoms
          l_end = (pjx != pix) ? 0 : 1
          // for all atom(s) in 'pix'
          for ( aid = 0 ; aid < pix_i->patoms - l_end ; aid ++ )
          {
            // set pointer to the first atom to test
            at_i = & pix_i->pix_atom[aid];
            lstart = (pjx != pix) ? 0 : aid + 1
            // for all atom(s) in 'pix_j'
            for ( bid = lstart ; bid < pix_j->patoms ; bid ++ )
            {
              // set pointer to the second atom to test
              at_j = & pix_j->pix_atom[bid];
              // evaluate interatomic distance
              Dij = evaluate_distance (use_pbc, at_i, at_j);
              if ( Dij.length < cutoff_squared )
              {
                // this is a first neighbor bond !
              }
            }
          }
        }
      }
      // store that pixel 'pix' was tested
      pix_i->tested = TRUE;
    }
  }
}
