/*
 * Copyright 2006-2008 The MZmine Development Team
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

package net.sf.mzmine.data.impl;

import java.text.Format;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Vector;

import net.sf.mzmine.data.CompoundIdentity;
import net.sf.mzmine.data.Peak;
import net.sf.mzmine.data.PeakListRow;
import net.sf.mzmine.io.RawDataFile;
import net.sf.mzmine.main.MZmineCore;
import net.sf.mzmine.project.MZmineProject;

/**
 * 
 */
public class SimplePeakListRow implements PeakListRow {

    private Hashtable<String, Peak> peaks;
    private Hashtable<String, Peak> originalPeaks;
    private HashSet<CompoundIdentity> identities;
    private CompoundIdentity preferredIdentity;
    private String comment;
    private int myID;
    private float maxDataPointIntensity = 0;

    public SimplePeakListRow(int myID) {
        this.myID = myID;
        peaks = new Hashtable<String, Peak>();
        originalPeaks = new Hashtable<String, Peak>();
        identities = new HashSet<CompoundIdentity>();
        preferredIdentity = CompoundIdentity.UNKNOWN_IDENTITY;
    }

    /**
     * @see net.sf.mzmine.data.PeakListRow#getID()
     */
    public int getID() {
        return myID;
    }

    /*
     * Return peaks assigned to this row
     */
    public Peak[] getPeaks() {
        return peaks.values().toArray(new Peak[0]);
    }

    /*
     * Return opened raw data files with a peak on this row
     */
    public RawDataFile[] getRawDataFiles() {
        MZmineProject project=MZmineCore.getCurrentProject();
        Vector <RawDataFile>rawDataFilesInPeak= new Vector<RawDataFile>();
        for (String fileName:peaks.keySet()){
        	rawDataFilesInPeak.add(project.getDataFile(fileName));
        }
        return rawDataFilesInPeak.toArray(new RawDataFile[0]);
    }

    /*
     * Returns peak for given raw data file
     */
    public Peak getPeak(RawDataFile rawData) {
        return peaks.get(rawData.getFileName());
    }

    public void setPeak(RawDataFile rawData, Peak p) {
        peaks.put(rawData.getFileName(), p);
    }

    /**
     * @see net.sf.mzmine.data.PeakListRow#getOriginalPeakListEntry(net.sf.mzmine.io.RawDataFile)
     */
    public Peak getOriginalPeakListEntry(RawDataFile rawData) {
        return originalPeaks.get(rawData.getFileName());
    }

    public void addPeak(RawDataFile rawData, Peak original, Peak current) {
        peaks.put(rawData.getFileName(), current);
        originalPeaks.put(rawData.getFileName(), original);
        if (current.getRawDataPointMaxIntensity() > maxDataPointIntensity)
            maxDataPointIntensity = current.getRawDataPointMaxIntensity();
    }

    /*
     * Returns average normalized M/Z for peaks on this row
     */
    public float getAverageMZ() {
        float mzSum = 0.0f;
        Enumeration<Peak> peakEnum = peaks.elements();
        while (peakEnum.hasMoreElements()) {
            Peak p = peakEnum.nextElement();
            mzSum += p.getMZ();
        }
        return mzSum / peaks.size();
    }

    /*
     * Returns average normalized RT for peaks on this row
     */
    public float getAverageRT() {
        float rtSum = 0.0f;
        Enumeration<Peak> peakEnum = peaks.elements();
        while (peakEnum.hasMoreElements()) {
            Peak p = peakEnum.nextElement();
            rtSum += p.getRT();
        }
        return rtSum / peaks.size();
    }

    /*
     * Returns number of peaks assigned to this row
     */
    public int getNumberOfPeaks() {
        return peaks.size();
    }

    public String toString() {
        StringBuffer buf = new StringBuffer();
        Format mzFormat = MZmineCore.getDesktop().getMZFormat();
        Format timeFormat = MZmineCore.getDesktop().getRTFormat();
        buf.append("#" + myID + " ");
        buf.append(mzFormat.format(getAverageMZ()));
        buf.append(" m/z @");
        buf.append(timeFormat.format(getAverageRT()));
        if (preferredIdentity != null)
            buf.append(" " + preferredIdentity.getCompoundName());
        if ((comment != null) && (comment.length() > 0))
            buf.append(" (" + comment + ")");
        return buf.toString();
    }

    /**
     * @see net.sf.mzmine.data.PeakListRow#getComment()
     */
    public String getComment() {
        return comment;
    }

    /**
     * @see net.sf.mzmine.data.PeakListRow#setComment(java.lang.String)
     */
    public void setComment(String comment) {
        this.comment = comment;
    }

    /**
     * @see net.sf.mzmine.data.PeakListRow#addCompoundIdentity(net.sf.mzmine.data.CompoundIdentity)
     */
    public void addCompoundIdentity(CompoundIdentity identity) {
        identities.add(identity);
        if (preferredIdentity == null)
            setPreferredCompoundIdentity(identity);
    }

    /**
     * @see net.sf.mzmine.data.PeakListRow#getCompoundIdentities()
     */
    public CompoundIdentity[] getCompoundIdentities() {
        return identities.toArray(new CompoundIdentity[0]);
    }

    /**
     * @see net.sf.mzmine.data.PeakListRow#getPreferredCompoundIdentity()
     */
    public CompoundIdentity getPreferredCompoundIdentity() {
        return preferredIdentity;
    }

    /**
     * @see net.sf.mzmine.data.PeakListRow#setPreferredCompoundIdentity(net.sf.mzmine.data.CompoundIdentity)
     */
    public void setPreferredCompoundIdentity(CompoundIdentity identity) {
        preferredIdentity = identity;
    }

    /**
     * @see net.sf.mzmine.data.PeakListRow#getDataPointMaxIntensity()
     */
    public float getDataPointMaxIntensity() {
        return maxDataPointIntensity;
    }

    public boolean hasPeak(Peak peak) {
        return peaks.containsValue(peak);
    }

}
