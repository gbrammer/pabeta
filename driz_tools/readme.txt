These scripts are intended to provide automation for much of the drizzle
reduction process.  The scripts are to be run in a specific order.

Each target should have it's own directory.  Each directory should contain
the FLT (or FLC if available) images for the target.  The scripts should be
run in the following order:

1. vis_driz.py
2. reg_drzs.py
3. back_wcs.py
4. calc_bounds.py
5. fin_driz.py

The scripts do the following:
1. Make visit level drizzle products for each filter. (f*dr?.fits)
2. Use visit level drizzle products to determine alignment shifts.
3. Propagates the shifts from the visit level products to individual exposures
4. Calculates bounding box for all images in RA/Dec space (so all final products
    have same dimensions).
5. Make final, full combined mosaics with desired scales etc. (F*dr?.fits)

Some of the scripts work off of various command line arguments call them with
the -h flag for more information.
