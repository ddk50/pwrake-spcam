# Rake Workflow for SDFRED1 (Subaru SprimeCam Data Reduction)

# SDFRED1 Reference:
# - Yagi et al., 2002, AJ, 123, 66
# - Ouchi et al., 2004, ApJ, 611, 660

# Requires SExtractor (2.5.0)
# - Bertin & Arnouts 1996, A&AS, 117, 393


dir1=Dir.glob(Dir.pwd+"/sextractor-*/src").last
dir2=Dir.glob(Dir.pwd+"/sdfred*/bin").last

ENV['PATH'] = "#{dir1}:#{dir2}:"+ENV['PATH']

require "rake/clean"
require "./sdfred.rb"

SDFREDSH=Dir.glob("sdfred*/sdfredSH").last

SFITS = FileList["SUPA*.fits"]



# Step 1: Basic Data Inspection and Renaming of Data Frames
#
# $ namechange.csh [raw fits file list]
#
#   raw fits file list = names of raw data files
#
# Prior to data reduction, it is useful to rename the data files using
# information associated with sensible parameters e.g., date of
# observations, exposure, and component CCD.
#
# The filename such as SUPA... is changed as
# H[Date][type][ID]_[chipname].fits.
#
# where the date, YYMMDD, is one day prior to DATE-OBS, and
# corresponds to Hawaiian Standard Time (HST) of the first half of the
# night.
#
# ID is the frame serial number of the observation day of each
# type(bias,dark,object). It should be noted that both target
# object(s) and standard star(s) have the same ID of "object".
#
# Example:
#
# $ cd object/
#
# enters the directory of object frames
#
# $ ls -1 SUPA*.fits > namechange.lis
#
# Notice that the option of "ls -1" is "minus one".
#
# $ namechange.csh namechange.lis
#
# The namechange.lis should be like that
#
# $ cat namechange.lis
# SUPA00204610.fits
# SUPA00204611.fits
# SUPA00204612.fits
# ...
#
# Result
# The files are renamed as follows;
#
# H030425object025_si001s.fits
# H030425object025_si002s.fits
# H030425object025_si005s.fits
# ...
#
# If you make symbolic links to the files in ../spcam_training_data/,
# they still points to SUPA***, but their names show up as
# H030425object025_si001s.fits. That is OK.
#
# Note that each CCD of Suprime-Cam has a name:
#
#     -------- AG probe location ----------
#     w67c1 w6c1 si006s si002s w7c3
#     w93c2 w9c2 si005s si001s w4c5

require "date"

HFITS_HASH = {}
HFITS = []

SFITS.each do |src|
  n = `getkey S_UFNAME #{src}`
  d = Date.parse(`getkey DATE-OBS #{src}`) - 1
  #name = "H#{d.year%100}#{d.mon}#{d.day}#{n}"
  name = "H%02d%02d%02d%s" % [d.year%100, d.mon, d.day, n]
  # puts "rename #{src} -> #{name}"
  HFITS_HASH[name] = src
  # File.rename( src, name )
  HFITS << name
end

HFITS.sort!



# Step 2: Subtraction overscan and bias
#
#  $ overscansub.csh [overscansub.lis]
#
#    overscansub.lis = list of raw data files
#
#  The script overscansub.csh issues a command that subtracts the
# median value of the overscan region in each line, and trims the
# overscan region off from the frame. The Suprime-Cam CCDs typically
# have an overscan level of about 10000 ADU.
#
#  Since the CCDs in Suprime-Cam have very little bias pattern, our
# experience suggests that subtracting overscan should suffice for
# many cases.
#
# Example:
#
# $ ls -1 H*.fits > ovserscansub.lis
# $ overscansub.csh ovserscansub.lis
#
# The overscansub.lis should be like that
#
# $ cat overscansub.lis
# H030425object025_si001s.fits
# H030425object025_si002s.fits
# H030425object025_si005s.fits
# ...
#
# Checkpoints:
#
#  * Compare the values of original frame and overscan subtracted
#    image. The latter should be about 10000 ADU smaller
#  * Checking the file headers with a procedure like IRAF's imhead
#    (cl> imhead H*.fits) should show that the image sizes are smaller
#    after removing the overscan region.

#HFITS = FileList["H[0-9]*.fits"]

TO_RHFITS = []
HFITS.each do |image|
  target = "To_R#{image}"
  file target => HFITS_HASH[image] do |t|
    sh "osmed4 -dtype=FITSFLOAT -bzero=0 -bscale=1 -pixignr=-32768 #{t.prerequisites[0]} #{t.name}"
  end
  TO_RHFITS << target
end

task :overscansub => TO_RHFITS
task :step2 => :overscansub

CLEAN.include(TO_RHFITS)



# Step 3: Making flat field frames
#
#  $ mask_mkflat_HA.csh [mkflat.lis]  [base name]  [lower value]  [upper value]
#
#    mkflat.lis = list of files to use to make flats
#    base name = basename for the flats
#    lower value = minimum value to accept (0.5 is recommended)
#    upper value = maximum value to accept (1.3 is recommended)
#
# The script mask_mkflat_HA.csh creates a flat from files with
# objects. The flat file is used to correct the difference in
# sensitivities between pixels in a frame. Areas vignetted by the AG
# probe is masked out, normalized, and a median of them is taken.
#
# There are three basic types of flats: sky flats (blank fields),
# twilight flats, and dome flats. Sky flats usually give the best
# result. This is because the optical path of the light for the night
# sky as well as the wavelength dependency are similar to those in
# your scientific frames. In fact, the target frames can be used to
# produce flats as long as there are no large objects that extend
# several hundred pixels in the frames.
#
# Example:
#
# $ ls -1 To_RH*.fits >mkflat.lis
# $ mask_mkflat_HA.csh mkflat.lis obj 0.5 1.3
#
# This is an example of making a sky flat by combining the object frames
#
# After running the script, there should be flat files for the 10 CCDs.
#
# obj_mflat_si001s.fits
# obj_mflat_si002s.fits
# ...
#
# The mflat files should have values around unity and should have a
# smooth pattern without much local structure. U-band and bands redder
# than z will have more structure than other bands. However, any local
# variations should be continuous. If there are abrupt changes in the
# flat values, consider creating another flat after eliminating
# (possible) bad exposures. Note that the value -32768 is the blank
# value used in SDFRED, and it is OK.
#
# Note 1:
#
# In principle, a flat can be produced with a minimum of three
# exposures. However, the smaller the number of frames used, the
# larger the noise and residual effects of objects in the frame. We
# recommend using at least six frames, ideally over 20 frames, to
# produce a (sensible) flat --- especially if you attempt to make a
# flat from sky frames. In the sample data, the flat is created from
# five exposures, yielding noticeable (but acceptable for robust
# reduction) error.
#
# Note 2:
#
# Keep in mind that the users should not mix the different types of
# flat exposures to make a flat frame because the background
# illuminations have intrinsically different slopes.  For example,
# SDFRED may produce flats with discontinuous stripes when applied to
# frames with different illumination patterns. This is due to the
# algorithm used in SDFRED.
#
# Note 3:
#
# The SDFRED command uses a parameter file to mask out known bad
# columns and hot pixels. The default parameter file is optimized for
# data taken after August 2002 (Messia V). To reduce data taken before
# this, follow the procedures below.
#
# $ cd sdfred20080610
#
# For data taken before March 2001 (old CCD):
#
# $ cp sdfredSH/mask_mkflat_HA/blankmap_oldCCD/* sdfredSH/mask_mkflat_HA/blankmap/
#
# For data taken between April 2001 and August 2002 (Messia III):
#
# $ cp sdfredSH/mask_mkflat_HA/blankmap_messiaIII/* sdfredSH/mask_mkflat_HA/blankmap/
#
# To return to the default:
#
# $ cp sdfredSH/mask_mkflat_HA/blankmap_messiaV/* sdfredSHmask_mkflat_HA/blankmap/
#
# Note 4:
#
# In the sample data, the flat is created in object/ directory, and
# will be shared for reducing standard-star data (by copying it into
# the standard/ directory). In many cases, several targets taken with
# an identical filter configuration can be mixed to create flats. If
# this is the case, it would be prudent to create another directory
# ../flat/. Subsequently, make symbolic links of To_R* files in all
# targets, and execute any commands in flat/ directory. (Note that the
# standard star frames that have shorter exposure times must not be
# included in the list.)
#
#
# e.g.)
# ---(work directory root) - object1/
#                          - object2/
#                          - object3/
#                          - flat/
#                          - standard/

MASK_MKFLAT_HA_SRC="#{SDFREDSH}/mask_mkflat_HA/mask_mkflat_HA.sex"
MASK_MKFLAT_HA_CONF="tmp_mask_mkflat_HA.sex"

file( MASK_MKFLAT_HA_CONF => MASK_MKFLAT_HA_SRC ) do
  sh "sed '/USER_BACK_SIZE/ d' #{MASK_MKFLAT_HA_SRC} > #{MASK_MKFLAT_HA_CONF}"
end

CLEAN.include MASK_MKFLAT_HA_CONF

# mask_mkflat_HA.csh mkflat.lis obj 0.5 1.3
HEAD  = "obj"
LOWER = 0.5
UPPER = 1.3
BLANK = -32768


HTO_RHFITS = []
TO_RHFITS.each do |image|
  mask_mkflat_sex = "mask_mkflat_HA_" + image.sub(/\.fits$/, ".sex")
  check_mask_mkflat_fits = "check_mask_mkflat" + image
  mask_mkflat_cat = "mask_mkflat" + image.sub(/\.fits$/, ".cat")

  #CLEAN.include mask_mkflat_sex
  CLEAN.include check_mask_mkflat_fits
  CLEAN.include mask_mkflat_cat

  #file( mask_mkflat_sex => image ) {|t|
  file( check_mask_mkflat_fits => [image, MASK_MKFLAT_HA_CONF] ) {|t|
    sex_mesh_s = 128 #determine the mesh size for SExtractor
    xpix = `getkey NAXIS1 #{image}`.to_i
    ypix = `getkey NAXIS2 #{image}`.to_i
    x_residual = xpix % sex_mesh_s
    y_residual = ypix % sex_mesh_s
    x_multiple = xpix / sex_mesh_s
    y_multiple = ypix / sex_mesh_s

    #for x_axis
    if x_residual==0
      x_add = 0
    else
      if x_residual <= sex_mesh_s/2.0
        x_add = x_residual / x_multiple + 1
      else
        x_add = (x_residual - sex_mesh_s) / x_multiple + 1
      end
    end
    sex_mesh_x = x_add + sex_mesh_s

    #for y_axis
    if y_residual == 0
      y_add = 0
    else
      if y_residual <= sex_mesh_s/2.0
        y_add = y_residual / y_multiple + 1
      else
        y_add = (y_residual - sex_mesh_s) / y_multiple + 1
      end
    end
    sex_mesh_y = y_add + sex_mesh_s

    printf( "#sexmesh %s %d %d %d %d    %.2f %.2f\n",
            image, sex_mesh_x, sex_mesh_y,
            xpix%sex_mesh_x, ypix%sex_mesh_y,
            1-(xpix%sex_mesh_x)/sex_mesh_x,
            1-(ypix%sex_mesh_y)/sex_mesh_y )

    opt  = "-c #{MASK_MKFLAT_HA_CONF} "
    opt += "-CATALOG_NAME /dev/null "
    opt += "-BACK_SIZE '#{sex_mesh_x},#{sex_mesh_y}' "
    opt += "-CHECKIMAGE_NAME #{check_mask_mkflat_fits} "
    sh "sex #{t.prerequisites[0]} #{opt} > /dev/null"

  }
  CLEAN.include check_mask_mkflat_fits

  tmp_hole_fits = "tmp_hole_" + image
  file tmp_hole_fits => [check_mask_mkflat_fits] do |t|
    sh "uppercut -imin=-10 -imax=0.0 -pixignr=-32768 #{t.prerequisites[0]} #{t.name}"
  end
  CLEAN.include tmp_hole_fits

  himage = "h#{image}"
  file himage => [image, tmp_hole_fits] do |t|
    sh "arithimg #{t.prerequisites[0]} - #{t.prerequisites[1]} #{t.name}"
  end
  HTO_RHFITS << himage

end

CLEAN.include HTO_RHFITS


task :mask => HTO_RHFITS
task :hto => HTO_RHFITS


AHTO_RHFITS = []
HTO_RHFITS.each do |himage|
  # shimasku note 5/17, 2002
  #(1) si006s, si002s, w6c1
  #  mask the region :   j > 60*AGX - 2300
  #(2) w67c1, w7c3
  #  mask the region :   j > 60*AGX - 2000

  ahimage = "a#{himage}"

  file ahimage => himage do |t|
    agx = `getkey S_AG-X #{t.prerequisites[0]}`.to_f
    #p agx
    j_limit =
    case himage
    when /si006s|si002s|w6c1/
      60 * agx - 2300
    when /w67c1|w7c3/
      60 * agx - 2000
    else
      5000
    end
    sh "mask_for_AGX #{t.prerequisites[0]} #{t.name} #{j_limit} #{BLANK}"
    #rm h${image}
  end
  AHTO_RHFITS << ahimage

end
CLEAN.include AHTO_RHFITS

task :ahto => AHTO_RHFITS


NAHTO_RHFITS = []
AHTO_RHFITS.each do |ahimage|

  nahimage = "n#{ahimage}"

  file nahimage => ahimage do |t|
    #appricable for data taken before and after 2001 April (for old and new CCDs)
    wmd_col1 = 30
    wmd_col2 = 1850
    wmd_row1 = 30
    wmd_row2_upperCCD = 3400
    wmd_row2_lowerCCD = 4060

    if /si002s|si006s|w6c1|w67c1|w7c3/ =~ t.name
      wmd_row2 = wmd_row2_upperCCD
    else
      wmd_row2 = wmd_row2_lowerCCD
    end

    sh "sync"
    sh "wmediandiv2 #{t.prerequisites[0]} #{wmd_col1} #{wmd_col2} #{wmd_row1} #{wmd_row2} #{t.name}"
    #rm ah${image}
  end

  NAHTO_RHFITS << nahimage

end
CLEAN.include NAHTO_RHFITS

task :nahto => NAHTO_RHFITS


MNAHTO_RHFITS = []
NAHTO_RHFITS.each do |nahimage|
  mnahimage = "m#{nahimage}"
  file mnahimage => nahimage do |t|
    sh "uppercut -imin=#{LOWER} -imax=#{UPPER} -pixignr=-32768 #{t.prerequisites[0]} #{t.name}"
    #rm nah${image}
  end
  MNAHTO_RHFITS << mnahimage
end
CLEAN.include MNAHTO_RHFITS

task :mnahto => MNAHTO_RHFITS



CHIPS = %w[si001s si002s si005s si006s w67c1 w6c1 w93c2 w9c2 w4c5 w7c3]

TMP_MKFLAT_FITS = []
TMP_MKFLAT_LIST = []
MFLAT = []

CHIPS.each do |chip|
  prereq = MNAHTO_RHFITS.select{|x| x=~/#{chip}/}
  blank = "tmp_mkflat_#{chip}.fits"

  file blank => prereq do |t|
    rejection_sigma = ''
    number_of_rejection = ''
    list = t.name+".list"
    File.open(list,"w"){|f| t.prerequisites.each{|x| f.puts x}}
    sh "sync"
    sh "mcomb2 #{list} #{t.name} #{rejection_sigma} #{number_of_rejection}"
  end

  TMP_MKFLAT_FITS << blank
  TMP_MKFLAT_LIST << blank+".list"

  blankmap = "#{SDFREDSH}/mask_mkflat_HA/blankmap/blankmap_spcamred_mflat_#{chip}"
  mflat = "#{HEAD}_mflat_#{chip}.fits"

  file mflat => [blank, blankmap] do |t|
    sh "blank2 #{t.prerequisites.join ' '} 60000 -32768 #{t.name}"
    if !File.readable?(t.name)
      copy t.prerequisites[0], t.name
    end
  end

  MFLAT << mflat
end

CLEAN.include MFLAT
CLEAN.include TMP_MKFLAT_FITS
CLEAN.include TMP_MKFLAT_LIST

task :mkflat => MFLAT
task :step3 => :mkflat



# Step 4: Flat Fielding
#
#  $ ffield.csh [ffiled_mf.lis]  [ffield_im.lis]
#
#    ffield_mf.lis = list of flats to be used
#    ffield_im.lis = list of (overscan subtracted) files to be flat fielded
#
# This command corrects the pixel-to-pixel variation in sensitivity,
# and the effect of vignetting of the telescope optics.
#
# Example:
#
# $ ls -1 obj_mflat*.fits > ffield_mf.lis
# $ ls -1 To_*.fits > ffield_im.lis
# $ ffield.csh ffield_mf.lis ffield_im.lis
#
# After the flat fielding, the background in each file should be
# almost flat. The circular illumination pattern seen in the raw data
# at the edge of the focal plane (w67c1, w93c2, w7c3, w4c5) should be
# nearly gone at this stage. If the flat is poor, you will see a
# low-level (several percent of variation) illumination
# pattern. However, even if you are not happy with the residual
# pattern, you can safely ignore it and proceed to Step 5 as long as
# you don't wish to achieve very high accuracy in photometry (i.e.,
# down to several percent). This is because such an illumination
# pattern will be subtracted at the later, at Step 8.
#
# On the other hand, if your scientific goals require a rather high
# accuracy in photometry --- down to several percent, consider going
# back to Step 3 to make another flat. The "ffield.csh" used in Step 4
# will help you to see the difference in the patterns appeared in the
# frames listed in mkflat.lis. We suggest eliminating the frames that
# have different patterns from those of the object frames from
# mkflat.lis in the Step 3.

FTORHFITS = []
TO_RHFITS.each do |image|
  chip = CHIPS.find{|x| image.include? x}
  mflat = MFLAT.find{|x| x.include? chip}
  fimage = "f"+image
  file fimage => [image, mflat] do |t|
    sh "arithimg #{t.prerequisites[0]} / #{t.prerequisites[1]} #{t.name}"
  end
  FTORHFITS << fimage
end

CLEAN.include FTORHFITS

task :flatfield => FTORHFITS
task :step4 => :flatfield



# Step 5: Distortion correction and atmospheric dispersion correction
#
#  $ distcorr.csh [distcorr.lis]
#
#    distcorr.lis = list of (flat fielded) files to be corrected
#
# The script distcorr.csh corrects the field distortion due to the
# telescope optics, and the differential atmospheric dispersion. The
# input frames are assumed to be flat-fielded images. The corrections
# are based on the airmass and other values recorded in the FITS
# header.
#
# Example:
#
# $ ls -1 fTo_RH030*.fits >distcorr.lis
# $ distcorr.csh distcorr.lis
#
# After the distortion correction, diagonal patterns may show up in
# the background. This is due to the fact that fractional pixel shifts
# smooth out the noise while integer pixel shifts leave the original
# noise characteristics intact.


GFTORHFITS = []
FTORHFITS.each do |image|
  gimage = "g#{image}"
  file gimage => [image] do |t|
    sh "sync"
    sh "distcorr5.sh #{t.prerequisites[0]} #{t.name}"
  end
  GFTORHFITS << gimage
end

CLEAN.include GFTORHFITS

task :distortion => GFTORHFITS
task :step5 => GFTORHFITS



# Step 6: Measurement of PSF size
#
#  $ fwhmpsf_batch.csh [fwhmpsf_batch.lis] [max number of objects]
#    [min peak flux] [max peak flux] [min FWHM] [max FWHM]
#
#    fwhmpsf_batch.lis = list of images to check PSF
#    max number of objects = the number of stars to use to measure
#                            the PSF in each image
#    min peak flux = minimum peak flux of stars to use
#    max peak flux = maximum peak flux of stars to use
#    min FWHM = minimum FWHM of stars to use
#    max FWHM = maximum FWHM of stars to use
#
# Before coadding, equalization of the PSF is required. The script
# fwhmpsf_batch.csh is used to determine an appropriate target PSF for
# a list of images. The script measures the FWHM of the PSF in several
# images. The script outputs a log and a histogram (exposure by
# exposure) to the standard output.
#
# Example:
#
# $ ls -1 gfTo_RH03042*.fits > fwhmpsf_batch.lis
# $ fwhmpsf_batch.csh fwhmpsf_batch.lis 50 2000 30000 2.0 7.0 \
#  > fwhmpsf_batch.log
#
# The log file should contain entries such as:
#
# gfTo_RH030425object025_si001s.fits   3.60   1 6 15 13 0
# gfTo_RH030425object025_si002s.fits   3.80   1 1 20 16 0
# gfTo_RH030425object025_si005s.fits   3.60   2 13 19 0 0
# ...
# 3.3 |
# 3.4 |
# 3.5 |**
# 3.6 |**
# 3.7 |*
# 3.8 |**
# 3.9 |**
# 4.0 |
# 4.1 |*
#
# The log format of each line is as follows:
#
#       [name of image]
#
#       [mean FWHM of PSF]
#
#       [number of objects with within 0.1'' of mean PSF-0.2]
#
#       [# within 0.1'' of mean PSF-0.1]
#
#       [# within 0.1'' of mean PSF]
#
#       [# within 0.1'' of mean PSF+0.1]
#
#       [# within 0.1'' of mean PSF+0.2]
#
# An ASCII histogram following the log illustrates the distribution of mean PSFs.
#
# Some training is required to determine the target PSF size.
#
# Here, the sample images have FWHM values centered around 3.7
# pixels. Adopting Target=4.1 would be the most conservative selection
# if you are interested in measuring fluxes of objects, but not in
# their structure. Target=3.9 would be a choice, if you believe the
# difference between 4.1 and 3.9 is negligible. Alternatively, if you
# have plenty of data, you may want to adopt a more strict threshold
# to exclude images having degraded PSFs (e.g, all data whose PSF are
# larger than 3.7).  The selection depends on your scientific goals.

# Note 1:
#
#  $ fwhmpsf.csh [image file] [max number of objects]
#    [min peak flux] [max peak flux] [min FWHM] [max FWHM]
#
#    image file = the image to check PSF
#    max number of objects = the number of stars to use to measure
#                            the PSF in each image
#    min peak flux = minimum peak flux of stars to use
#    max peak flux = maximum peak flux of stars to use
#    min FWHM = minimum FWHM of stars to use
#    max FWHM = maximum FWHM of stars to use
#
# Use the script fwhmpsf.csh to find the PSF of a single image. The
# parameters are the same as for fwhmpsf_batch.csh. Just supply the
# name of a image rather than a list of images.
#
# Example:
# $ fwhmpsf.csh gfTo_RH030425object025_si001s.fits 50 2000 30000 2.0 7.0
#
# This produces output like:
#
# gfTo_RH030425object025_si001s.fits   3.60   1 6 15 13 0
#
# This output indicates that the image
# gfTo_RH030425object025_si001s.fits has a PSF FWHM of 3.6 pixels.

# Note 2:
#
#  $ starselect.csh [image name][max number of objects][min peak flux]
#      [max peak flux][min FWHM][max FWHM][output file]
#
#    image name = name of image to check
#    max number of objects = the number of stars to use to measure the PSF
#    min peak flux = minimum peak flux of stars to use
#    max peak flux = maximum peak flux of stars to use
#    min FWHM = minimum FWHM of stars to use
#    max FWHM = maximum FWHM of stars to use
#    output file = name of file with location of selected stars
#
# The script starselect.csh is useful for searching for the
# appropriate parameters ([max number of objects] [min peak flux] [max
# peak flux] [min FWHM] and [max FWHM]) for selecting stellar objects
# in an image.
#
# $ starselect.csh gfTo_RH030425object025_si001s.fits 50 2000 \
#   30000 2.0 7.0 output.reg
#
# The script will produce an output file (output.reg) that contains
# the location of stellar objects satisfying the given criteria. The
# output is formatted so that the stellar objects are plotted with
# green circles when you plot using saoimage-ds9. If the majority of
# the selected objects are "real stellar objects" (stars for many
# cases), then the parameters are appropriate for psf_match for a
# given image. If you realize that the quality of the data varies
# image by image, determine whether or not a single set of parameters
# can be applied for whole the data set.  If it cannot, it is better
# to run psfmatch_batch multiple times using the appropriate criteria
# for each subset of data.
#
# Using saoimage-ds9 is the easiest way to display an image and
# overlay the location of the selected stars.
#
# $ ds9 gfTo_RH030425object025_si001s.fits
#
# Select "Region", "Load", and select output.reg. Then green circles
# will be overlaid on the image. If more than half of the objects
# selected are stellar objects, the parameters you adopted are
# appropriate

# Note 3:
#
# The scripts fwhmpsf_batch.psf, starselect.csh and
# psfmatch_batch.csh(next step), may not work in crowded fields. In
# such fields, it may be necessary to estimate the PSF manually. If
# this is the case, to obtain the same results as the psfmatch script,
# each image that has a PSF more than 0.1'' smaller than the target
# PSF should be Gaussian smoothed with a Gaussian that has a
# sigma=sqrt(PSF_target^2 - PSF_image^2)/2.35482.  The psfmatch
# program iterates around the estimated value of the sigma until the
# results converges to the best matched final PSF.


# $ ls -1 gfTo_RH03042*.fits > fwhmpsf_batch.lis
# $ fwhmpsf_batch.csh fwhmpsf_batch.lis 50 2000 30000 2.0 7.0 > fwhmpsf_batch.log

FWHMPSF_BATCH_LIS = GFTORHFITS.select{|x| /gfTo_RH03042/=~x }

MF_number_of_obj = 50
MF_min_flux = 2000
MF_max_flux = 30000
MF_min_fwhm = 2.0
MF_max_fwhm = 6.0

FWHMPSF_CAT = []
FWHMPSF_CHK = []

FWHMPSF_BATCH_LIS.each do |item|
  catlg = File.basename(item,".fits")+".fwhm.cat"

  file catlg => item do|t|
    # fwhmpsf.csh
    image = t.prerequisites[0]
    chkimg = "check_" + image
    sh "sex #{image} -c #{SDFREDSH}/fwhmpsf/fwhmpsf.sex -CATALOG_NAME #{t.name} -CHECKIMAGE_NAME #{chkimg}> /dev/null"
  end

  FWHMPSF_CAT << catlg
  FWHMPSF_CHK << "check_"+item
end

CLEAN.include FWHMPSF_CAT
CLEAN.include FWHMPSF_CHK


file "tmp_fwhmpsf_batch1" => FWHMPSF_CAT do |t|

  tmp1 = t.prerequisites.map do |cat|
    image = File.basename(cat,".cat")+".fits"
    tmp2 = select_obj( cat, MF_min_flux, MF_max_flux,
                       MF_min_fwhm, MF_max_fwhm, MF_number_of_obj )
    fwhmpsf = histmax( tmp2, MF_min_fwhm, MF_max_fwhm, 0.1 )
    upper_fwhmpsf = fwhmpsf+0.2
    lower_fwhmpsf = fwhmpsf-0.2
    d = mkhist( tmp2, lower_fwhmpsf, upper_fwhmpsf, 0.1 )
    [image, fwhmpsf, d[0][1], d[1][1]]
  end

  open("tmp_fwhmpsf_batch1","w") do |f|
    tmp1.each do |x|
      f.puts x.join(' ')
    end
  end
end

file "tmp_fwhmpsf_batch2" => "tmp_fwhmpsf_batch1" do |t|
  tmp = []
  open(t.prerequisites[0]) do |f|
    while s=f.gets
      a = s.split
      tmp << a[1].to_f if a[1]
    end
  end

  open(t.name, "w") do |f|
    target_fwhm = make_a_histogram(tmp, 0.1, f)
    f.printf "target_fwhm = %.1f\n", target_fwhm
    puts "********************"
    printf " target_fwhm = %.1f\n", target_fwhm
    puts "********************"
  end
end

CLEAN.include "tmp_fwhmpsf_batch1"
CLEAN.include "tmp_fwhmpsf_batch2"

task :step6 => "tmp_fwhmpsf_batch2"



# Step 7: Equalize the PSF size
#
#  $ psfmatch_batch.csh [psfmatch_batch.lis] [max number of objects]
#      [min peak flux] [max peak flux] [min FWHM] [max FWHM] [target FWHM]
#
#    psfmatch_batch.lis = the list of images to match to a single PSF
#    max number of objects = the number of stars to use to measure the
#                            PSF in each image
#    min peak flux = minimum peak flux of stars to use
#    max peak flux = maximum peak flux of stars to use
#    min FWHM = minimum FWHM of stars to use
#    max FWHM = maximum FWHM of stars to use
#    target FWHM = FWHM to smooth all the data to
#
# The script psfmatch_batch.csh attempts to match the PSF of all
# images to be combined to a predetermined target FWHM. Images with
# PSFs smaller than the target (within a small range) are Gaussian
# smoothed, other images are simply copied. The target PSF should
# represent the typical PSF for the exposure, having the worst (i.e.,
# the largest) PSF among the exposures to be combined.
#
# The command prints a log to the standard output with the following columns:
#
# [name of psf_matched_image]
# [FWHM of PSF after matching]
# [number of objects with within 0.1'' of target PSF-0.2]
# [# within 0.1'' of  target PSF-0.1]
# [# within 0.1'' of target PSF]
# [# within 0.1'' of  target PSF+0.1]
# [# within 0.1'' of  target PSF+0.2]
#
# The log can be used to check whether or not the PSF matching worked
# properly. To determine whether or not the command has ended
# successfully, make sure that the number of objects falling into (1)
# the bin of +/-0.1" of the final FWHM must be the peak of the
# distribution, and (2) outside bins e.g., PSF+0.2 arcsecond and
# PSF-0.2 are significantly small.
#
#   *
#   *
#   *
#   *
#  **
#  **
#  ***
#  ***
#  ***
# *****
#
# If this is not the case, you are strongly encouraged to check the PSF manually.
#
# Example:
#
# $ ls -1 gfTo_RH03042*.fits > psfmatch_batch.lis
# $ psfmatch_batch.csh psfmatch_batch.lis 50 2000 30000 2.0 7.0 3.7 > psfmatch_batch.log &
#
# The command produces output like:
#
# pgfTo_RH030429object017_si001s.fits   3.70   0 5 38 6 0
#
# The IRAF task "imexam" is handy for checking PSFs. (Display image;
# cl> imexam image.fits; place cursor above a star; type "r" or "a" to
# measure FWHM)
#
# Keep in mind that each software routine may be using different
# fitting algorithms and may return different FWHM values. SDFRED
# adopts FWHM values generated by SExtractor, which are different from
# those produced by IRAF's imexam task. The purpose of checking with
# IRAF is not to find an exact match in the FWHM values, but to
# confirm that the output images have comparable PSF sizes after the
# matching.
#
# The appropriate parameter values for psfmatch will change depending
# on the quality of the data. Different bandpasses, integration times,
# and weather conditions will require different parameters.


# $ ls -1 gfTo_RH03042*.fits > psfmatch_batch.lis
# $ psfmatch_batch.csh psfmatch_batch.lis 50 2000 30000 2.0 7.0 3.7 > psfmatch_batch.log &

PSFMATCH_BATCH_LIS = GFTORHFITS.select{|x| /gfTo_RH03042/=~x }

#psfmatch_batch.csh psfmatch_batch.lis 50 2000 30000 2.0 7.0 3.7
P_number_of_obj = 50
P_min_flux = 2000
P_max_flux = 30000
P_min_fwhm = 2.0
P_max_fwhm = 7.0
# P_target_fwhm = 3.7

#load "Rakefile.psfmatch"

PSFMATCH_CAT = []
PSFMATCH_OUT = []
PSFMATCH_LOG = []

PSFMATCH_BATCH_LIS.each do |image_in|
  image_cat = File.basename(image_in,".fits")+".cat"
  PSFMATCH_CAT << image_cat

  image_out = "p"+image_in
  logfile = File.basename(image_in,".fits")+".log"

  file image_out => [image_in,"tmp_fwhmpsf_batch2"]  do |t|
    img = t.prerequisites[0]
    out = t.name # "p"+img
    target_fwhm = read_target_fwhm(t.prerequisites[1])
    target_fwhm = 3.7
    puts "target_fwhm=#{target_fwhm}"
    log = File.basename(img,".fits")+".log"
    sh "ruby psfmatch.rb #{img} #{P_number_of_obj} #{P_min_flux} #{P_max_flux} #{P_min_fwhm} #{P_max_fwhm} #{target_fwhm} #{out} #{log}"
  end

  file logfile => image_out

  PSFMATCH_OUT << image_out
  PSFMATCH_LOG << logfile
end

CLEAN.include PSFMATCH_CAT
CLEAN.include PSFMATCH_OUT
CLEAN.include PSFMATCH_LOG


file "psfmatch_batch.log" => PSFMATCH_LOG do |t|
  tmp_psfmatch_batch2 = []
  t.prerequisites.each do |fn|
    open(fn) do |f|
      while s=f.gets
        if /^pg/=~s
          a = s.split
          tmp_psfmatch_batch2 << a[1].to_f if a.size > 1
        end
      end
    end
  end

  open(t.name, "w") do |f|
    make_a_histogram(tmp_psfmatch_batch2, 0.1, f)
  end

end

CLEAN.include "psfmatch_batch.log"

task :step7 => "psfmatch_batch.log"


def read_target_fwhm(file)
  open(file) do |f|
    while s=f.gets
      return $1.to_f if /target_fwhm\s*=\s*([\d.]+)/ =~ s
    end
  end
end


# Step 8: Subtracting the Sky Background
#
#  $ skysb.csh [skysb.lis] [sky-mesh]
#
#    skysb.lis = list of images to sky subtract
#    sky-mesh = size of mesh for determining sky values
#
# The script skysb.csh (1) computes a mesh pattern that represents the
# sky background, (2) interpolates the pattern, and (3) subtracts it
# from the image. The script creates a grid --- referred to as
# "sky-mesh size squares" --- on the image with a grid spacing having
# the half of the "sky-mesh" size. An appropriate sky-mesh size will
# be selected for each mesh, and assigned to the pixel located at the
# center of the mesh. After rejecting the outliers, the sky values for
# other pixels will be given by interpolating bilineary from the
# surrounding meshes. Note that the sky-mesh size must be selected at
# least twice the largest object in interest due to the Nyquist
# sampling theorem.
#
# Example:
# $ ls -1 pgfTo_RH03042*.fits > skysb.lis
# $ skysb.csh skysb.lis 64 > skysb.log
#
# Once the sky background level is subtracted, the background in an
# image should be around zero without a spatial gradient. If there is
# an extended object(s) spreading over a large fraction of the image,
# the algorithm will most likely fail. Subtraction of sky background
# in crowded fields requires special data handling and you will need
# to estimate the sky background manually.

SKYSB_LIS = PSFMATCH_OUT.select{|x| /pgfTo_RH03042/=~x }
SKYSB_OUT = []
SKYSB_MESH = 64

SKYSB_LIS.each do |img_in|
  img_out = "s"+img_in
  file img_out => img_in do |t|
    blankvalue = -32768
    sh "skysb3b -imin=-1000 -imax=32500 -pixignr=#{blankvalue} -mesh=#{SKYSB_MESH} #{t.prerequisites[0]} #{t.name}"
  end
  SKYSB_OUT << img_out
end

CLEAN.include SKYSB_OUT
task :step8 => SKYSB_OUT



# Step 9: Masking the AG Shade
#
#  $ mask_AGX.csh [mask_AGX.lis]
#
#    mask_AGX.lis = list of files to mask
#
# The script mask_AGX.csh will mask areas vignetted by the AG
# (Auto-Guider) probe by the value -32768. The script should only
# affect the top few hundred rows of the data from chips w671, w6c1,
# si006s, si002s, and w7c3. Other files are not affected.
#
# Example:
# $ ls -1 spgfTo_RH03042*.fits > mask_AGX.lis
# $ mask_AGX.csh mask_AGX.lis
#
# Although only half the CCDs are potentially affected by the AG
# probe, the input file list should include all the object files so
# that files with the same naming convention exist to make list-making
# for subsequent steps easier.

MASK_AGX_LIS = SKYSB_OUT # .select{|x| /pgfTo_RH03042/=~x }
MASK_AGX_OUT = []

MASK_AGX_LIS.each do |img_in|
  img_out = "A"+img_in
  file img_out => img_in do |t|
    image = t.prerequisites[0]
    blank = -32768
    agx = `getkey S_AG-X #{image}`.to_i
    case image
    when /si006s|si002s|w6c1/
      j_limit = agx * 60 - 2300
    when /w67c1|w7c3/
      j_limit = agx * 60 - 2000
    else
      j_limit = 5000
    end
    sh "mask_for_AGX #{image} #{t.name} #{j_limit} #{blank}"
  end
  MASK_AGX_OUT << img_out
end

CLEAN.include MASK_AGX_OUT
task :step9 => MASK_AGX_OUT



# Step 10: Masking Bad Pixels
#
# Data in some pixels may be corrupted due to instrument trouble
# and/or other problems which may have occurred during the
# observation. Such regions must be common among the exposures (i.e.,
# they are not time variable), and should be masked accordingly. For
# instance, we suggest masking the background areas where flattening
# fails and systematically deviates from zero. If plenty of exposures
# cover the observed region, we suggest not spending much time with
# this step. This is because outliers will be rejected automatically
# in Step 12.
#
# The SDFRED1 package offers three methods --- linear, circular, and
# rectangular regions --- to specify regions to be masked for
# eliminating bad pixels. Here, "Linear region" connects the two
# points (x1,y1) - (x2,y2), extends the line to the edges of the
# image, and masks the pixels within "width" from the line. The
# "circular region" masks the pixels in a circle. The "rectangular
# regions" masks rectangular regions aligned to the pixel coordinate.
#
# Linear region
#
#  $ line_bank [input image] [x1] [y1] [x2] [y2] [width]
#     [blank value] [output image]
#
#    input image = name of image to mask
#    x1 = x coordinate of start of line
#    y1 = y coordinate of start of line
#    x2 = x coordinate of end of line
#    y2 = y coordinate of end of line
#    width = width of line
#    blank value = mask value (usually -32768)
#    output image = name of masked image
#
# The script line_bank masks a linear structure such as satellite trails.
#
# Example:
#
# $ line_blank AspgfTo_RH030425object025_si001s.fits \
#   88 112 1940 837 30 -32768 lAspgfTo_RH030425object025_si001s.fits
#
# The example masks line which crosses (88,112) and (1940,837) and around 30 pixels width.
# Circular region
#
#  $ circular_blanks [input image] [blanklist] [blank value] [output image]
#
#    input image = name of image to mask
#    blank list = a text file describing the x and y coordinates, as
#                 well as the radii of the areas to be masked
#    blank value = mask value (usually -32768)
#    output image = name of masked image
#
# The script circular_blanks masks circular regions.
#
# Example:
#
# $ circular_blanks lAspgfTo_RH030425object025_si001.fits \
#   blanklist -32768 clAspgfTo_RH030425object025_si001s.fits
#
# where blanklist looks like:
#
# $ cat blanklist
# 365 1835 80
# 1202 3582 100
#
# The two lines correspond to circle of (x,y,r)=(365,1835,80) and (x,y,r)=(1202,3582,100).
# Rectangular regions
#
#  $ blank.csh [blank list]
#
#    blank list = list of images to be masked
#
# For each image, xxx.fits, in the blank list the script blank.csh
# will look for a file named blankmap_xxx in the same directory, and
# mask rectangular regions specified in the file to -32768. Each line
# in the file blank_xxx should contain the x and y coordinates of two
# opposite corners of a rectangular area.
#
# The IRAF routine imexam is useful for getting the coordinates. (cl>
# imexam; press "b" at two corners to define a rectangle; the
# coordinates of the corners will be printed to the screen in the
# order of x1 x2 y1 y2.)
#
# Example:
#
# $ ls -1 AspgfTo_RH03042*.fits > blank.lis
# $ blank.csh blank.lis
#
# Mask files have been included for a subset of images,
#
#       blankmap_AspgfTo_RH030425object025_si002s
#
#       ...
#
# These files have entries like:
#
# $ cat blankmap_AspgfTo_RH030425object025_si002s
# 1974 2034 2356 2634
# 1528 1804 4024 4070
# ...
#
# The script masks the regions specified in the corresponding
# blankmap_xxx file. If the blankmap_xxx file does not exist, the
# script will simply copy the image file to the output.
#
# Note that the blankmap_* parameter files included in the sample
# dataset are from the Subaru Deep Field project and are not
# necessarily suitable for the sample dataset. If all of the sample
# files are applied to the sample data, areas that shouldn't normally
# be masked will be masked. These files are strictly for practice
# use. You may wish to create your own mask files to apply to the
# sample data.


BLANK_LIS = MASK_AGX_OUT
BLANK_OUT = []

BLANK_LIS.each do |img_in|
  img_out = "b"+img_in
  blankmap = "blankmap_" + File.basename(img_in,".fits")
  if File.file?(blankmap)
    file img_out => [img_in,blankmap] do |t|
      sh "sync; blank2 #{t.prerequisites.join(' ')} 60000 -32768 #{t.name}"
    end
  else
    file img_out => img_in do |t|
      sh "sync; cp #{t.prerequisites[0]} #{t.name}"
    end
  end
  BLANK_OUT << img_out
end

CLEAN.include BLANK_OUT
task :step10 => BLANK_OUT



# Step 11: Estimating Alignment and Scaling
#
#  $ makemos.csh [makemos.lis] [starsel nskysigma] [starsel npix]
#                [starsel peakmin] [starsel peakmax]
#                [aperture phot radius in pix] [output mos-file name]
#
#    makemos.lis = list of images to align
#    starsel nskysigma = signal to noise ratio of objects per pixel
#                        to use for alignment
#    starsel npix = number of continuous pixels with [starsel nskysigma]
#                   to identify object
#    starsel peakmin = minimum value of peak pixel of alignment stars
#    starsel peakmax = maximum value of peak pixel of alignment stars
#    aperture phot radius in pix = radius to use for aperture photometry
#    output mos-file name = file to record alignment and scaling
#
# Signal-to-noise ratio (S/N) can be improved by combining multiple
# images (if you have them) to produce a final image. The script
# makemos.csh determines the shifts, rotations, and flux scales of
# different images. The script identifies stellar objects in each
# image and determines the the shifts, rotations, and flux scale from
# stellar objects common to multiple images. The first image in the
# list is used as the reference image.
#
# Example:
#
# $ ls -1 bAspgfTo_RH03042*.fits > makemos.lis
# $ makemos.csh makemos.lis 5 20 500 10000 10 all.mos > makemos.log
#
# The script will print to the standard output the number of stellar
# objects selected for alignment and scaling.
#
#      ...
#      selected stars = 721
#      ...
#
# The script is likely to fail if the number of selected stars per
# image is either small (< 30) or very large (>1000). Optimizing key
# parameters such as [starsel nskysigma], [starsel npix], [starsel
# peakmin], and [starsel peakmax] will help the script to select
# appropriate stellar objects.
#
# The best parameters for selecting objects in this step may be
# different from PSF measurement for many cases. This is because a
# different underlying algorithm is employed to find a wider range of
# objects to determine relative positions and flux scaling that works
# over a range of fluxes.
#
# Example:
#
# $ cat all.mos
# bAspgfTo_RH030425object025_si001s.fits 0.000000 0.000000 0.000000 1.000000
# bAspgfTo_RH030425object025_si002s.fits 1.601005 4088.551744 -0.000275 0.989406
# bAspgfTo_RH030425object025_si005s.fits -2118.787371 1.104977 -0.000282 0.983195
# ...
#
# As shown in the above, you will see five parameters (i.e., columns)
# in the output *.mos file: the name of the image, the x offset, the y
# offset, the counter clockwise rotation (radian), and flux ratio. If
# each result has four output parameters followed by the image name,
# the alignment or/and scaling have successfully finished. If the
# alignment and scaling have failed, the output file may not be
# produced at all, miss some parameters, or have unreasonable values.
#
# We -- the SDFRED support team -- have been making continuous efforts
# to provide users more sophisticated method(s) that examines
# all.mos. Realizing the situation, however, we wish to share the
# following tips:
#
#    1. Inspecting the final image created in the next step must be
#    done. However, bear in mind that it is not the ultimate
#    method. If the number of exposure is large, it is difficult to
#    detect some small defects by visual inspection in the final
#    image.
#    2. It is always a good idea to make a plot of the 2nd vs. 3rd
#    columns stored in all.mos. The result shows the relative position
#    of each shot, and represents the dither pattern as well as the
#    chip positions. If there is a large leap in value, the matching
#    has failed.
#    3. The distances between CCD chips should be almost constant. (A
#    slight difference may exist due to atmospheric dispersion between
#    chips.) If the distances between any arbitrarily chosen chip
#    pairs for the same exposure (e.g., between si001s and w67c1) has
#    changed significantly exposure by exposure, the data of the
#    corresponding exposure would be incorrect.
#    4. The fifth column of *.mos (relative flux) of a chip should be
#    almost proportional to the exposure time, if sky condition is
#    photometric. (It is affected by atmospheric extinction (airmass),
#    however)
#
# In the next step, each image is converted with the data in all.mos as follows;
#
#
# x_mos =   cos(theta) x  -sin(theta) y + x_local
# y_mos =   sin(theta) x  +cos(theta) y + y_local
#
# Note 1:
# If you don't need to combine, you can skip Steps 11 and 12, and end
# the reduction. If you intend to combine images toward more than two
# fields, make sure that these data have been taken contiguously. If
# this is not the case, Step 11 will fail.
#

# $ ls -1 bAspgfTo_RH03042*.fits > makemos.lis
# $ makemos.csh makemos.lis 5 20 500 10000 10 all.mos > makemos.log

MOS_nskysigma = 5
MOS_npix = 20
MOS_peakmin = 500
MOS_peakmax = 10000
MOS_apphotrad = 10
MOS_outmos = "all.mos"

MES_LIS = []

MAKEMOS_LIS = BLANK_OUT

MAKEMOS_LIS.each do |image|
  mes = image + ".mes"
  file mes => image do |t|
    sh "sync"
    sh "starsel2 -nskysigma=#{MOS_nskysigma} -npix=#{MOS_npix} -aratiomin=0.4 -peakmin=#{MOS_peakmin} -peakmax=#{MOS_peakmax} -aprad=#{MOS_apphotrad} -outmes=#{t.name} #{t.prerequisites[0]} > /dev/null"
  end
  MES_LIS << mes
end
CLEAN.include MES_LIS

# shotmatch7a.sh 30 `cat shotmatch.lis` > log_shotmatch

if ! File.file? "overlap_lis"
  OVERLAP_TASK_LIST = []
  OVERLAP_LIST = []
  count = 0
  MES_LIS.each do |ames|
    if /To_R(H\d{6}.*)\.mes$/ =~ ames
      aorg = HFITS_HASH[$1]
    end
    MES_LIS.each do |bmes|
      if /To_R(H\d{6}.*)\.mes$/ =~ bmes
        borg = HFITS_HASH[$1]
      end
      if ames != bmes
        OVERLAP_TASK_LIST << n = "match#{count}"
        task( n => [aorg,borg] ) do |t|
          sh "overlap2 #{t.prerequisites.join(' ')}" do |res,status|
            #if status.to_i == 0
            if res
              OVERLAP_LIST << x = "#{t.name} #{ames} #{bmes}"
              puts x
            end
          end
        end
        count += 1
      end
    end
  end

  if ! Rake.application.top_level_tasks.include?("clean")

    file( "overlap_lis" => OVERLAP_TASK_LIST ) do |t|
      open(t.name,"w") do |f|
        OVERLAP_LIST.each{|x| f.puts x}
      end
    end.invoke

  end
end


CLEAN.include "overlap_lis"

SM_nmin = 4 # for Suprime
SM_nmax = 30
#MATCHSINGLE="match_single5"

MATRIX_LIS=[]
open("overlap_lis") do |f|
  while s = f.gets
    mdat, ames, bmes = s.split
    mdat += '.dat'

    file( mdat => [ames,bmes] ) do |t|
      a, b = t.prerequisites
      shotmatch( a, b, SM_nmin, SM_nmax, t.name )
    end

    MATRIX_LIS << mdat
  end
end

CLEAN.include MATRIX_LIS


file "matrix.dat" => MATRIX_LIS do |t|
  open(t.name,"w") do |w|
    t.prerequisites.each do |fn|
      s = IO.read(fn).chomp
      w.puts s if /^bAspgf/ =~ s
    end
  end
end

CLEAN.include "matrix.dat"

MATCHSTACK="match_stack5"

file "all.mos" => "matrix.dat" do |t|
  sh "#{MATCHSTACK} #{t.prerequisites[0]} > #{t.name}.tmp"
  sh "cat #{t.name}.tmp | sort > #{t.name}"
end

CLEAN.include ["all.mos", "all.mos.tmp"]

task :step11 => "all.mos"



# Step 12: Combining
#
#  $ imcio2a [parameters] [mos file] [result image]
#
#  parameters = parameters that define the combining algorithm usually
#   "-dist_clip -nline=20 -dtype=FITSFLOAT -pixignr=-32768"
#  mos file = file containing the alignment and scaling values
#             (output from makemos.csh)
#  result image= the name of the final image
#
# imcio2a combines the images into a final combined image using the
# output from makemos.csh (*.mos). Using the parameter -dist_clip will
# combine the images using a clipped mean algorithm.
#
# Example:
#
# $ imcio2a -dist_clip -nline=20 -dtype=FITSFLOAT -pixignr=-32768 all.mos all.fits
#
# The parameter -dist_clip can be replaced by -dist_med to get a
# weighted median combined image or -dist_add to use a weighted mean
# pixel values.
#
# Note that the header of the output image is incomplete. Use the
# first file listed in makemos.lis in the previous step as a reference
# header.
#
# Here are the meanings of the typical parameters:
#
#       -dist_clip : use a clipped mean algorithm for combining
#
#       -nline=20 : set the y direction buffer width to 20
#
#       -dtype=FITSFLOAT : make the output data floating point
#
#       -pixignr=-32768 : ignore pixels valued -32768
#
# For details and other optional parameters of imcio2a can be printed by
#
#  $ imcio2a -h
#

# $ imcio2a -dist_clip -nline=20 -dtype=FITSFLOAT -pixignr=-32768 all.mos all.fits


file "all.fits" => ["all.mos"]+BLANK_OUT do |t|
  if false # defined? Pwrake::GfarmSSH and conn = Thread.current[:connection]
    gfpwd = Pwrake::GfarmSSH.gf_pwd
    hostname = conn.host
    list = t.prerequisites.map{|x| "#{gfpwd}/#{x}"}.join(" ")
    cmd = "gfrep -m -N 1 -D #{hostname} #{list}"
    puts cmd
    system cmd
  end
  sh "imcio2a -dist_clip -nline=20 -dtype=FITSFLOAT -pixignr=-32768 #{t.prerequisites[0]} #{t.name}"
end

task :default => "all.fits"


# Reduction of Standard Object
#
# Steps 1S through 4S describe a typical procedure for reducing
# standard stars data. Since the physics behind this is the same as
# for reducing target objects, you can essentially repeat the
# procedure. Don't forget to work in the standard/ directory. The flat
# frames must be the same as those used for the objects, therefore be
# sure to copy them from the object/ directory.

# Step 1S Renaming
#
# $ namechange.csh [raw fits file list]
#
#   raw fits file list = names of raw data files
#
# Renaming is done in the standard/ directory
#
# Example:
#
# $ cd standard/
#
# enters into the directory of standard frames
#
# $ ls -1 SUPA*.fits > namechange.lis
# $ namechange.csh namechange.lis
#
# The namechange.lis should be like that
#
# $ cat namechange.lis
# SUPA00195120.fits
# SUPA00195121.fits
# SUPA00195122.fits
# ...
#
# Result
# The files are renamed as follows;
#
# H030330object044_si001s.fits
# H030330object044_si002s.fits
# H030330object044_si005s.fits
# ...

# Step 2S Overscan and bias subtraction
#
#  $ overscansub.csh [overscansub.lis]
#
#    overscansub.lis = list of raw data files
#
# In standard/ directory, overscan is subtracted from all the data as follows,
#
# Example:
#
# $ ls -1 H*.fits > overscansub.lis
# $ overscansub.csh overscansub.lis
#
# $ cat overscansub.lis
# H030330object044_si001s.fits
# H030330object044_si002s.fits
# H030330object044_si005s.fits
# ...
#
# and, To_RH030330object044_si001s.fits ... are created.

# Step 3S Flat fielding
#
#  $ ffield.csh [ffiled_mf.lis]  [ffield_im.lis]
#
#    ffield_mf.lis = list of flats to be used
#    ffield_im.lis = list of (overscan subtracted) files to be flat fielded
#
# The flat frames used in this step (ffield_mf.lis) must be identical
# to those used for the target(s) in order to cancel out the
# uncertainty in the normalization.
#
# Example:
#
# $ cp ../object/obj_mflat*.fits .
# $ ls -1 obj_mflat*.fits > ffield_mf.lis
# $ ls -1 To_RH*.fits > ffield_im.lis
# $ ffield.csh ffield_mf.lis ffield_im.lis
#
# and fTo_RH030330object044_si001s.fits ... are created.

# Step 4S Distortion correction and atmospheric dispersion correction
#
#  $ distcorr.csh [distcorr.lis]
#
#    distcorr.lis = list of (flat fielded) files to be corrected
#
# The distortion correction is required since it slightly changes the
# sizes of the pixels, yielding slightly different flux value(s).
#
# Example:
#
# $ ls -1 fTo_RH*.fits >distcorr.lis
# $ distcorr.csh distcorr.lis

# Step 5S Correction of relative flux scale among chips
#
# Recall that the relative flux "scale" between different CCD chips
# has not yet been corrected even after the Step 4S. For example, the
# sensitivity of w67c1 is about one half in 2002-2008/06 data. If
# there is a star that has 10000 ADU in w67c1, it should have ~ 20000
# ADU if it was observed with the other chip. The relative flux scale
# should be corrected according to the *.mos created in step (11) of
# the target object.
#
# This step is unnecessary if the standard star is only in
# si001s. Since standard stars are distributed in several chips in the
# sample data, this process is needed.
#
# In this step, the data is divided by a typical relative value of the
# chip to the reference chip. Currently (SDFRED ver1.*), the script
# for this step is not provided. Users should do this process
# manually.
#
# Acknowledgment
# The SDFRED team thanks Takehiko Wada (JAXA), Chiaki Ihara (JAXA),
# Hitoshi Hanami (Iwate Univ), Myungkook James Jee (JHU), Kazuaki Ota
# (NAOJ), Yasunori Sato (NAOJ), Ryosuke Yamauchi (Tohoku Univ.),
# Elinor Medezinski (Tel-Aviv Univ.), Dovi Poznanski (Tel-Aviv Univ.),
# Ben Cain (MIT), Tomoki Saito (Ehime Univ.), Sakurako Okamoto
# (Univ. of Tokyo), Alice Shapley (UCLA) and Naoki Yasuda (Univ. of
# Tokyo) for their comments and bug reports.
