import shutil
import os
import time

import numpy

from azcam.console import azcam
import azcam_testers
from .basetester import Tester


class DetCal(Tester):
    """
    Detector calibration routines to:
     - find and set video offsets
     - find exposure levels in DN and electrons at specified wavelengths
     - find system gains
     - read diode flux calibration data
    """

    def __init__(self):

        super().__init__("detcal")

        # offsets
        self.offsets = []
        self.old_offsets = []
        self.offset_filename = ""
        self.dsp_compiler = ""
        self.video_scale = -0.5  # ARC-48
        # self.video_scale=-3.3       # ARC-45
        self.bias_goal = 1000
        self.mean_count_goal = 10000
        self.zero_image = "test.fits"
        self.data_file = "detcal.txt"

        # calibration
        self.exposure_type = "flat"
        self.overwrite = 0  # True to overwrite old data
        self.wavelength_delay = 2  # seconds to delay after changing wavelengths
        self.zero_mean = []
        self.system_gain = []

        self.wavelengths = []  # list of list of wavelengths to calibrate
        self.exposure_times = {}  # list of dictionaries of {wavelength:initial guess et}
        self.mean_counts = {}  # list of dictionaries of {wavelength:Counts/Sec}
        self.mean_electrons = {}  # list of dictionaries of {wavelength:Electrons/Sec}

    def calibrate(self):
        """
        Use gain data to find offsets and gain.
        Take images at each wavelength to get count levels.
        If no wavelength are specified, only calibrate current wavelength
        """

        azcam.log("Running detector calibration sequence")

        # save pars to be changed
        impars = {}
        azcam.api.save_imagepars(impars)

        # create new subfolder
        if self.overwrite:
            if os.path.exists("detcal"):
                shutil.rmtree("detcal")
        startingfolder, subfolder = azcam.utils.make_file_folder("detcal")
        azcam.api.set_par("imagefolder", subfolder)
        azcam.utils.curdir(subfolder)

        azcam.api.set_par("imageincludesequencenumber", 1)  # don't use sequence numbers
        azcam.api.set_par("imageautoname", 0)  # manually set name
        azcam.api.set_par("imagetest", 0)  # turn off TestImage
        azcam.api.set_par("imageoverwrite", 1)

        # get gain and ROI
        self.system_gain = azcam.api.gain.get_system_gain()
        self.roi = azcam.utils.get_image_roi()

        gain = azcam.api.gain

        self.system_gain = gain.system_gain
        self.zero_mean = gain.zero_mean

        # clear device
        azcam.api.tests()

        self.mean_counts = {}
        self.mean_electrons = {}

        wavelengths = self.wavelengths

        # get flat at each wavelength
        for wave in wavelengths:

            # set wavelength
            wave = int(wave)
            wave1 = azcam.api.get_wavelength()
            wave1 = int(wave1)
            if wave1 != wave:
                azcam.log(f"Setting wavelength to {wave} nm")
                azcam.api.set_wavelength(wave)
                time.sleep(self.wavelength_delay)
                wave1 = azcam.api.get_wavelength()
                wave1 = int(wave1)
            azcam.log(f"Current wavelength is {wave1} nm")

            # take flat
            doloop = 1
            try:
                et = self.exposure_times[wave]
            except Exception:
                et = 1.0
            while doloop:
                azcam.api.set_par("imagetype", self.exposure_type)
                azcam.log(f"Taking flat for {et:0.3f} seconds")
                flatfilename = azcam.api.get_image_filename()
                azcam.api.expose(et, self.exposure_type, "detcal flat")

                # get counts
                bin1 = int(azcam.fits.get_keyword(flatfilename, "CCDBIN1"))
                bin2 = int(azcam.fits.get_keyword(flatfilename, "CCDBIN2"))
                binning = bin1 * bin2
                flatmean = numpy.array(azcam.fits.mean(flatfilename)) - numpy.array(self.zero_mean)
                flatmean = flatmean.mean()
                azcam.log("Mean signal at {wave} nm is {flatmean:0.0f} DN")

                if flatmean > self.mean_count_goal * 2.0:
                    et = et * (self.mean_count_goal / flatmean)
                    continue
                elif flatmean < self.mean_count_goal / 2.0:
                    et = et * (self.mean_count_goal / 2.0 / flatmean)
                    continue

                self.mean_counts[wave] = flatmean / et / binning
                self.mean_electrons[wave] = self.mean_counts[wave] * numpy.array(self.system_gain)

                self.mean_counts[wave] = self.mean_counts[wave].mean()
                self.mean_electrons[wave] = self.mean_electrons[wave].mean()
                doloop = 0

        # define dataset
        self.dataset = {
            "data_file": self.data_file,
            "wavelengths": self.wavelengths,
            "mean_electrons": self.mean_electrons,
            "mean_counts": self.mean_counts,
        }

        # write data file
        azcam.utils.curdir(startingfolder)
        self.write_datafile()

        self.valid = True

        # finish
        azcam.api.restore_imagepars(impars, startingfolder)
        azcam.log("detector calibration sequence finished")

        return

    def read_datafile(self, filename="default"):
        """
        Read data file and set object as valid.
        """

        super().read_datafile(filename)

        # convert types
        self.mean_counts = {int(k): v for k, v in self.mean_counts.items()}
        self.mean_electrons = {int(k): v for k, v in self.mean_electrons.items()}

        return
