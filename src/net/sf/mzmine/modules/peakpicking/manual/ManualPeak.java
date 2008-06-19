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

package net.sf.mzmine.modules.peakpicking.manual;

import java.util.TreeMap;

import net.sf.mzmine.data.DataPoint;
import net.sf.mzmine.data.Peak;
import net.sf.mzmine.data.PeakStatus;
import net.sf.mzmine.data.RawDataFile;
import net.sf.mzmine.util.CollectionUtils;
import net.sf.mzmine.util.MathUtils;
import net.sf.mzmine.util.PeakUtils;
import net.sf.mzmine.util.Range;

/**
 * This class represents a manually picked chromatographic peak.
 */
class ManualPeak implements Peak {

    private RawDataFile dataFile;

    // Raw M/Z, RT, Height and Area
    private float mz, rt, height, area;

    // Boundaries of the peak
    private Range rtRange, mzRange, intensityRange;

    private TreeMap<Integer, DataPoint> dataPointMap;

    /**
     * Initializes empty peak for adding data points
     */
    ManualPeak(RawDataFile dataFile) {
        this.dataFile = dataFile;
        dataPointMap = new TreeMap<Integer, DataPoint>();
    }

    /**
     * This peak is always a result of manual peak detection, therefore MANUAL
     */
    public PeakStatus getPeakStatus() {
        return PeakStatus.MANUAL;
    }

    /**
     * This method returns M/Z value of the peak
     */
    public float getMZ() {
        return mz;
    }

    /**
     * This method returns retention time of the peak
     */
    public float getRT() {
        return rt;
    }

    /**
     * This method returns the raw height of the peak
     */
    public float getHeight() {
        return height;
    }

    /**
     * This method returns the raw area of the peak
     */
    public float getArea() {
        return area;
    }

    /**
     * This method returns numbers of scans that contain this peak
     */
    public int[] getScanNumbers() {
        return CollectionUtils.toIntArray(dataPointMap.keySet());
    }

    /**
     * This method returns a representative datapoint of this peak in a given
     * scan
     */
    public DataPoint getDataPoint(int scanNumber) {
        return dataPointMap.get(scanNumber);
    }

    /**
     * This method returns all raw data points used to build this peak in a
     * given scan
     */
    public DataPoint[] getRawDataPoints(int scanNumber) {
        return new DataPoint[] { dataPointMap.get(scanNumber) };
    }

    public Range getRawDataPointsIntensityRange() {
        return intensityRange;
    }

    public Range getRawDataPointsMZRange() {
        return mzRange;
    }

    public Range getRawDataPointsRTRange() {
        return rtRange;
    }

    /**
     * @see net.sf.mzmine.data.Peak#getDataFile()
     */
    public RawDataFile getDataFile() {
        return dataFile;
    }

    public String toString() {
        return PeakUtils.peakToString(this);
    }

    /**
     * Adds a new data point to this peak
     * 
     * @param scanNumber
     * @param dataPoints
     * @param rawDataPoints
     */
    void addDatapoint(int scanNumber, DataPoint dataPoint) {

        float rt = dataFile.getScan(scanNumber).getRetentionTime();

        if (dataPointMap.isEmpty()) {
            rtRange = new Range(rt);
            mzRange = new Range(dataPoint.getMZ());
            intensityRange = new Range(dataPoint.getIntensity());
        } else {
            rtRange.extendRange(rt);
            mzRange.extendRange(dataPoint.getMZ());
            intensityRange.extendRange(dataPoint.getIntensity());
        }
        
        dataPointMap.put(scanNumber, dataPoint);

    }

    void finalizePeak() {

        // Trim the zero-intensity data points from the beginning and end
        while (!dataPointMap.isEmpty()) {
            int scanNumber = dataPointMap.firstKey();
            if (dataPointMap.get(scanNumber).getIntensity() > 0)
                break;
            dataPointMap.remove(scanNumber);
        }
        while (!dataPointMap.isEmpty()) {
            int scanNumber = dataPointMap.lastKey();
            if (dataPointMap.get(scanNumber).getIntensity() > 0)
                break;
            dataPointMap.remove(scanNumber);
        }

        // Check if we have any data points
        if (dataPointMap.isEmpty()) {
            throw (new IllegalStateException(
                    "Peak can not be finalized without any data points"));
        }

        // Get all scan numbers
        int allScanNumbers[] = CollectionUtils.toIntArray(dataPointMap.keySet());

        // Find the data point with top intensity and use its RT and height
        for (int i = 0; i < allScanNumbers.length; i++) {
            DataPoint dataPoint = dataPointMap.get(allScanNumbers[i]);
            float rt = dataFile.getScan(allScanNumbers[i]).getRetentionTime();
            if (dataPoint.getIntensity() > height) {
                height = dataPoint.getIntensity();
                this.rt = rt;
            }
        }

        // Calculate peak area
        area = 0f;
        for (int i = 1; i < allScanNumbers.length; i++) {

            // X axis interval length
            float previousRT = dataFile.getScan(allScanNumbers[i - 1]).getRetentionTime();
            float thisRT = dataFile.getScan(allScanNumbers[i]).getRetentionTime();
            float rtDifference = thisRT - previousRT;

            // Intensity at the beginning and end of the interval
            float previousIntensity = dataPointMap.get(allScanNumbers[i - 1]).getIntensity();
            float thisIntensity = dataPointMap.get(allScanNumbers[i]).getIntensity();
            float averageIntensity = (previousIntensity + thisIntensity) / 2;

            // Calculate area of the interval
            area += (rtDifference * averageIntensity);

        }

        // Calculate median MZ
        float mzArray[] = new float[allScanNumbers.length];
        for (int i = 0; i < allScanNumbers.length; i++) {
            mzArray[i] = dataPointMap.get(allScanNumbers[i]).getMZ();
        }
        this.mz = MathUtils.calcQuantile(mzArray, 0.5f);

    }

	public void setMZ(float mz) {
		// TODO Auto-generated method stub
		
	}
}
