/*
	Grid.c
		utility functions for the Vegas grid
		this file is part of Suave
		last modified 16 Jul 04 th
*/


static inline Grid *GetGrid()
{
  Grid *grid;

  if( vegasgridno_ < 1 || vegasgridno_ > MAXGRIDS ||
      (grid = gridptr_[vegasgridno_ - 1]) == NULL ||
      ((griddim_[vegasgridno_ - 1] != ndim_) ? free(grid), 1 : 0) ) {
    count bin, dim;

    Allocate(grid, ndim_*sizeof(Grid));

    for( bin = 0; bin < NBINS; ++bin )
      grid[0][bin] = (bin + 1)/(real)NBINS;
    for( dim = 1; dim < ndim_; ++dim )
      Copy(&grid[dim], &grid[0], 1);
  }

  return grid;
}

/*********************************************************************/

static inline PutGrid(Grid *grid)
{
  if( vegasgridno_ < 1 || vegasgridno_ > MAXGRIDS ) free(grid);
  else {
    gridptr_[vegasgridno_ - 1] = grid;
    griddim_[vegasgridno_ - 1] = ndim_;
  }
}

/*********************************************************************/

static void RefineGrid(Grid grid, real *margsum)
{
  real avgperbin, thisbin;
  Grid imp, newgrid;
  count bin, newbin;

  /* smooth the f^2 value stored for each bin */
  real prev = margsum[0];
  real cur = margsum[1];
  real norm = margsum[0] = .5*(prev + cur);
  for( bin = 1; bin < NBINS - 1; ++bin ) {
    creal s = prev + cur;
    prev = cur;
    cur = margsum[bin + 1];
    norm += margsum[bin] = (s + cur)/3.;
  }
  norm += margsum[NBINS - 1] = .5*(prev + cur);

  if( norm == 0 ) return;
  norm = 1/norm;

  /* compute the importance function for each bin */
  avgperbin = 0;
  for( bin = 0; bin < NBINS; ++bin ) {
    real impfun = 0;
    if( margsum[bin] > 0 ) {
      creal r = margsum[bin]*norm;
      avgperbin += impfun = pow((r - 1)/log(r), 1.5);
    }
    imp[bin] = impfun;
  }
  avgperbin /= NBINS;

  /* redefine the size of each bin */
  cur = 0;
  thisbin = 0;
  bin = -1;
  for( newbin = 0; newbin < NBINS - 1; ++newbin ) {
    while( thisbin < avgperbin ) {
      thisbin += imp[++bin];
      prev = cur;
      cur = grid[bin];
    }
    thisbin -= avgperbin;
    newgrid[newbin] = cur - (cur - prev)*thisbin/imp[bin];
  }
  Copy(grid, newgrid, NBINS - 1);
  grid[NBINS - 1] = 1;
}

/*********************************************************************/

static void Reweight(Grid *grid,
  creal *w, creal *f, creal *lastf, cCumulants *total)
{
  real margsum[NDIM][NBINS], scale[NCOMP];
  cbin_t *bin = (cbin_t *)lastf;
  count dim, comp;

  if( ncomp_ == 1 ) scale[0] = 1;
  else {
    for( comp = 0; comp < ncomp_; ++comp )
      scale[comp] = (total[comp].avg == 0) ? 0 : 1/total[comp].avg;
  }

  Zap(margsum);

  while( f < lastf ) {
    creal weight = *w++;
    for( comp = 0; comp < ncomp_; ++comp ) {
      creal fsq = Sq(*f++*weight*scale[comp]);
      if( fsq == 0 ) continue;
      for( dim = 0; dim < ndim_; ++dim )
        margsum[dim][bin[dim]] += fsq;
    }
    bin += ndim_;
  }

  for( dim = 0; dim < ndim_; ++dim )
    RefineGrid(grid[dim], margsum[dim]);
}

