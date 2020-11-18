import os

import azcam_testers

from azcam.console import azcam

from .basetester import Tester


class PocketPump(Tester):
    """
    Pocket pumping (trap) acquisition and analysis.
    """

    def __init__(self):

        super().__init__("pocketpump")

        self.exposure_type = "flat"
        self.exposure_time = 1.0
        self.exposure_level = -1
        self.number_cycles = 1000
        self.shift_rows = 1
        self.wavelength = 600.0

        self.data_file = "traps.txt"
        self.report_file = "traps"

        self.total_traps = 0  # should be -1

        self.means = []
        self.sdev = []
        self.temperatures = []
        self.pocketpump_timing_file = "STA3800C_config1.lod"  # example

    def acquire(self):
        """
        Make a pocket pump image sequence.
        """

        return self.acquire_archon()

    def acquire_arc(self):
        """
        Make a pocket pump image sequence.
        """

        azcam.log("Acquiring PocketPump sequence")

        # save pars to be changed
        impars = {}
        azcam.api.save_imagepars(impars)

        azcam.api.exposure.set_par(
            "imageroot", "pocketpump."
        )  # for automatic data analysis
        azcam.api.exposure.set_par(
            "imageincludesequencenumber", 1
        )  # use sequence numbers
        azcam.api.exposure.set_par("imageautoname", 0)  # manually set name
        azcam.api.exposure.set_par(
            "imageautoincrementsequencenumber", 1
        )  # inc sequence numbers
        azcam.api.exposure.set_par("imagetest", 0)

        # create and move to new subfolder
        currentfolder, newfolder = azcam.utils.make_file_folder("pocketpump")
        azcam.api.exposure.set_par("imagefolder", newfolder)

        # set wavelength so exposure time is correct
        azcam.log(f"Setting wavelength to {self.wavelength}")
        azcam.api.instrument.set_wavelength(self.wavelength)

        # set proper timing code, reset, and clear device
        timingfile_org = azcam.api.exposure.get_par("TimingFile")
        timcode = os.path.dirname(timingfile_org)
        timcode = os.path.join(timcode, self.pocketpump_timing_file)
        timcode = timcode.replace("\\", "/")
        azcam.log("PocketPump timing code will be: %s" % timcode)
        azcam.api.exposure.set_par("TimingFile", timcode)
        azcam.api.exposure.reset()
        azcam.api.exposure.test(2)

        # take a bias
        azcam.log("Taking first pocketpump bias")
        try:
            azcam.api.exposure.expose(0, "zero", "pocketpump first bias")
        except Exception as message:
            azcam.api.exposure.set_par("TimingFile", timingfile_org)
            azcam.api.exposure.reset()
            azcam.api.restore_imagepars(impars, currentfolder)
            raise message

        # setup exposure
        azcam.log("Pocketpump setup")
        azcam.api.exposure.begin_exposure(
            self.exposure_time, self.exposure_type, "pocket pump"
        )

        # integrate
        azcam.log("Pocketpump integration")
        azcam.api.exposure.integrate_exposure()

        # pocket pump
        azcam.log("Pocket pumping...")
        for _ in range(self.number_cycles):
            azcam.api.exposure.parshift(
                -1 * self.shift_rows
            )  # reverse par shift is a complete cycle

        # readout
        azcam.log("Readout")
        azcam.api.exposure.readout_exposure()

        # write image and display
        azcam.api.exposure.end_exposure()
        azcam.log("Pocket pumped exposure completed")

        # take a bias
        azcam.log("Taking second pocketpump bias")
        try:
            azcam.api.exposure.expose(0, "zero", "pocketpump second bias")
        except Exception as message:
            azcam.api.exposure.set_par("TimingFile", timingfile_org)
            azcam.api.exposure.reset()
            azcam.api.restore_imagepars(impars, currentfolder)
            return message

        # take a reference flat of same time
        azcam.log("Taking non-pumped reference exposure")
        try:
            azcam.api.exposure.expose(
                self.exposure_time, self.exposure_type, "pocket pump reference flat"
            )

        except Exception as message:
            azcam.api.exposure.set_par("TimingFile", timingfile_org)
            azcam.api.exposure.reset()
            azcam.api.restore_imagepars(impars, currentfolder)
            raise message

        # finish
        azcam.api.exposure.set_par("TimingFile", timingfile_org)
        azcam.api.exposure.reset()
        azcam.api.restore_imagepars(impars, currentfolder)
        azcam.log("PocketPump sequence finished")

        return

    def acquire_archon(self):
        """
        Make a pocket pump image sequence.
        """

        azcam.log("Acquiring PocketPump sequence")

        # save pars to be changed
        impars = {}
        azcam.api.save_imagepars(impars)

        azcam.api.exposure.set_par(
            "imageroot", "pocketpump."
        )  # for automatic data analysis
        azcam.api.exposure.set_par(
            "imageincludesequencenumber", 1
        )  # use sequence numbers
        azcam.api.exposure.set_par("imageautoname", 0)  # manually set name
        azcam.api.exposure.set_par(
            "imageautoincrementsequencenumber", 1
        )  # inc sequence numbers
        azcam.api.exposure.set_par("imagetest", 0)

        # create and move to new subfolder
        currentfolder, newfolder = azcam.utils.make_file_folder("pocketpump")
        azcam.api.exposure.set_par("imagefolder", newfolder)

        # set wavelength so exposure time is correct
        azcam.log(f"Setting wavelength to {self.wavelength}")
        azcam.api.instrument.set_wavelength(self.wavelength)

        # Get exposure time
        binning = 1
        if azcam.api.detcal.valid and self.exposure_level != -1:
            azcam.log("Using exposure_level")

            meanelectrons = azcam.api.detcal.mean_electrons

            et = self.exposure_level / meanelectrons[self.wavelength] / binning
        elif self.exposure_time != -1:
            azcam.log("Using ExposureTime")
            et = self.exposure_time
        else:
            raise azcam.AzcamError("could not determine exposure times")

        # flush well
        azcam.log("Flushing")
        azcam.api.exposure.test(2)

        # take a bias
        azcam.log("Taking first pocketpump bias")
        try:
            azcam.api.exposure.expose(0, "zero", "pocketpump first bias")
        except Exception as message:
            azcam.api.restore_imagepars(impars, currentfolder)
            raise message

        # turn on pocketpumping
        azcam.api.serverconn.rcommand("controller.set_pocket_pumping 1")

        # expose
        azcam.log("Taking pocketpump flat")
        azcam.api.exposure.expose(et, "flat", "pocketpump flat")

        # turn off pocketpumping
        azcam.api.serverconn.rcommand("controller.set_pocket_pumping 0")

        # take a bias
        azcam.log("Taking second pocketpump bias")
        try:
            azcam.api.exposure.expose(0, "zero", "pocketpump second bias")
        except Exception as message:
            azcam.api.restore_imagepars(impars, currentfolder)
            raise message

        # take a reference flat of same time
        azcam.log("Taking non-pumped reference exposure")
        try:
            azcam.api.exposure.expose(
                et, self.exposure_type, "pocket pump reference flat"
            )

        except Exception as message:
            azcam.api.restore_imagepars(impars, currentfolder)
            raise message

        # finish
        azcam.api.restore_imagepars(impars, currentfolder)
        azcam.log("PocketPump sequence finished")

        return

    def analyze(self):
        """
        Analys a pocket pump image sequence.
        Not yet supported.
        """

        self.total_traps = 0
        # self.grade='UNDEFINED'
        self.grade = "PASS"

        # define dataset
        self.dataset = {
            "data_file": self.data_file,
            "grade": self.grade,
            "total_traps": self.total_traps,
        }

        # write output files
        self.write_datafile()
        if self.create_reports:
            self.report()

        return
