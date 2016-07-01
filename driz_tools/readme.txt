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

IN ADDITION:
Two other scripts have been developed to increase overall data quality/usability:
1. ir_data_fixes.py
2. hsc_align.py

1. Flattens F110W ramp, flags extra bad pixels, and masks persistence
(Ideally run before step 1 in above, or have to redrizzle IR data)
2. Aligns final top level mosaics (F*dr?.fits) to Hubble Source Catalog, and updates input FLT headers accordingly.
(Run after step 5, or can just input HSC catalog into step 2)
Some of the scripts work off of various command line arguments call them with
the -h flag for more information.


Things that need to be done:
* Run ir_data_fixes on everything, and redrizzle accordingly
    If already astrometrized to HSC, need to run calc_bounds again,
    and then probably best to redrizzle the whole stack for a given target?
* Find good way to rename WCS of top level mosaics if flts already aligned to HSC
