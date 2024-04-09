_unused = 0x0001
no_data = 0x0002
bad_pix = 0x0004
ccd_gap = 0x0008
telluric = 0x0010

seg_has_badpixels = 0x0020

starlight_masked_pix = 0x0100
starlight_failed_run = 0x0200

d3d_screw = 0x1000

# Compound flags
no_obs = no_data | bad_pix | ccd_gap | telluric | d3d_screw | seg_has_badpixels
no_starlight = starlight_masked_pix | starlight_failed_run
