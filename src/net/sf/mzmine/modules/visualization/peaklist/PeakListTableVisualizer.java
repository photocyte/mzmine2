/*
 * Copyright 2006-2007 The MZmine Development Team
 * 
 * This file is part of MZmine.
 * 
 * MZmine is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 * 
 * MZmine is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * MZmine; if not, write to the Free Software Foundation, Inc., 51 Franklin St,
 * Fifth Floor, Boston, MA 02110-1301 USA
 */

package net.sf.mzmine.modules.visualization.peaklist;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.util.logging.Logger;

import net.sf.mzmine.data.ParameterSet;
import net.sf.mzmine.data.PeakList;
import net.sf.mzmine.main.MZmineCore;
import net.sf.mzmine.main.MZmineModule;
import net.sf.mzmine.userinterface.Desktop;
import net.sf.mzmine.userinterface.Desktop.MZmineMenu;
import net.sf.mzmine.util.logging.JCommonLogHandler;

import org.jfree.report.JFreeReportBoot;

public class PeakListTableVisualizer implements MZmineModule, ActionListener {

    private Desktop desktop;
    private PeakListTableParameters parameters;
    private static PeakListTableVisualizer myInstance;

    private Logger logger = Logger.getLogger(this.getClass().getName());

    /**
     * @see net.sf.mzmine.main.MZmineModule#initModule(net.sf.mzmine.main.MZmineCore)
     */
    public void initModule() {

        this.desktop = MZmineCore.getDesktop();
        myInstance = this;

        // boot the JFreeReport library and register our logging handler,
        // to get a rid of JCommon debug messages on the console
        JFreeReportBoot.getInstance().start();
        JCommonLogHandler.register();

        parameters = new PeakListTableParameters();

        desktop.addMenuItem(MZmineMenu.VISUALIZATION,
                "Aligned peak list table", this, null, KeyEvent.VK_A, false,
                true);

    }

    /**
     * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
     */
    public void actionPerformed(ActionEvent e) {

        PeakList peakLists[] = desktop.getSelectedPeakLists();

        if (peakLists.length == 0) {
            desktop.displayErrorMessage("Please select aligned peak list");
            return;
        }

        for (PeakList peakList : peakLists) {

            logger.finest("Showing a new aligned peak list table view");

            PeakListTableWindow alignmentResultView = new PeakListTableWindow(
                    peakList);
            desktop.addInternalFrame(alignmentResultView);

        }
    }

    /**
     * @see net.sf.mzmine.main.MZmineModule#toString()
     */
    public String toString() {
        return "Alignment result table visualizer";
    }

    /**
     * @see net.sf.mzmine.main.MZmineModule#getParameterSet()
     */
    public PeakListTableParameters getParameterSet() {
        return parameters;
    }

    /**
     * @see net.sf.mzmine.main.MZmineModule#setParameters(net.sf.mzmine.data.ParameterSet)
     */
    public void setParameters(ParameterSet parameterValues) {
        parameters = (PeakListTableParameters) parameterValues;
    }

    public static PeakListTableVisualizer getInstance() {
        return myInstance;
    }

}