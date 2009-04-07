/*
 * Copyright 2006-2009 The MZmine 2 Development Team
 *
 * This file is part of MZmine 2.
 *
 * MZmine 2 is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * MZmine 2 is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * MZmine 2; if not, write to the Free Software Foundation, Inc., 51 Franklin
 * St, Fifth Floor, Boston, MA 02110-1301 USA
 */
package net.sf.mzmine.project.impl;

import java.io.ByteArrayInputStream;



import java.io.ByteArrayOutputStream;
import org.jfree.xml.util.Base64;
import org.xml.sax.Attributes;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Hashtable;
import java.util.TreeMap;
import java.util.Vector;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;
import java.util.zip.ZipOutputStream;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;
import net.sf.mzmine.data.ChromatographicPeak;
import net.sf.mzmine.data.DataPoint;
import net.sf.mzmine.data.PeakIdentity;
import net.sf.mzmine.data.PeakList;
import net.sf.mzmine.data.PeakListAppliedMethod;
import net.sf.mzmine.data.PeakListRow;
import net.sf.mzmine.data.PeakStatus;
import net.sf.mzmine.data.RawDataFile;
import net.sf.mzmine.data.impl.SimpleChromatographicPeak;
import net.sf.mzmine.data.impl.SimpleDataPoint;
import net.sf.mzmine.data.impl.SimplePeakIdentity;
import net.sf.mzmine.data.impl.SimplePeakList;
import net.sf.mzmine.data.impl.SimplePeakListAppliedMethod;
import net.sf.mzmine.data.impl.SimplePeakListRow;
import net.sf.mzmine.data.impl.SimpleScan;
import net.sf.mzmine.main.mzmineclient.MZmineCore;
import net.sf.mzmine.modules.io.xmlexport.PeakListElementName;
import net.sf.mzmine.project.MZmineProject;
import net.sf.mzmine.util.Range;
import org.dom4j.Document;
import org.dom4j.DocumentFactory;
import org.dom4j.Element;
import org.dom4j.io.OutputFormat;
import org.dom4j.io.XMLWriter;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

public class PeakListSerializer extends DefaultHandler {

	private ZipOutputStream zipOutputStream;
	private ZipInputStream zipInputStream;
	public static DateFormat dateFormat = new SimpleDateFormat(
			"yyyy/MM/dd HH:mm:ss");
	private Hashtable<RawDataFile, Integer> dataFilesIDMap;
	private TreeMap<Integer, RawDataFile> buildingArrayRawDataFiles;
	private int parsedRows;
	private StringBuffer charBuffer;
	private Vector<String> appliedProcess;
	private PeakList buildingPeakList;
	private RawDataFileImpl buildingRawDataFile;
	private SimplePeakListRow buildingRow;
	private int peakColumnID,  rawDataFileID,  quantity;
	private double mass,  rt,  height,  area;
	private int[] scanNumbers;
	private double[] retentionTimes,  masses,  intensities;
	private String peakStatus,  peakListName,  name,  formula,  identificationMethod,  identityID;
	private boolean preferred;
	private String dateCreated;
	private Range rtRange,  mzRange;
	private boolean peakListFlag = false;
	private boolean scanFlag = false;
	private boolean mzPeakFlag = false;

	public PeakListSerializer(ZipOutputStream zipOutputStream) {
		this.zipOutputStream = zipOutputStream;
		dataFilesIDMap = new Hashtable<RawDataFile, Integer>();
	}

	public PeakListSerializer(ZipInputStream zipInputStream) {
		this.zipInputStream = zipInputStream;
	}

	public void savePeakList(PeakList peakList) throws IOException {
		Element newElement;
		Document document = DocumentFactory.getInstance().createDocument();
		Element saveRoot = document.addElement(PeakListElementName.PEAKLIST.getElementName());

		// <NAME>
		newElement = saveRoot.addElement(PeakListElementName.NAME.getElementName());
		newElement.addText(peakList.getName());

		// <PEAKLIST_DATE>
		String dateText = "";
		if (((SimplePeakList) peakList).getDateCreated() == null) {
			dateText = ((SimplePeakList) peakList).getDateCreated();
		} else {
			Date date = new Date();
			dateText = dateFormat.format(date);
		}
		newElement = saveRoot.addElement(PeakListElementName.PEAKLIST_DATE.getElementName());
		newElement.addText(dateText);

		// <QUANTITY>
		newElement = saveRoot.addElement(PeakListElementName.QUANTITY.getElementName());
		newElement.addText(String.valueOf(peakList.getNumberOfRows()));

		// <PROCESS>
		PeakListAppliedMethod[] processes = peakList.getAppliedMethods();
		for (PeakListAppliedMethod proc : processes) {
			newElement = saveRoot.addElement(PeakListElementName.PROCESS.getElementName());
			newElement.addText(proc.getDescription());
		}

		// <RAWFILE>
		RawDataFile[] dataFiles = peakList.getRawDataFiles();

		for (int i = 1; i <= dataFiles.length; i++) {
			newElement = saveRoot.addElement(PeakListElementName.RAWFILE.getElementName());
			newElement.addAttribute(
					PeakListElementName.ID.getElementName(), String.valueOf(i));
			// <NAME>
			newElement = saveRoot.addElement(PeakListElementName.NAME.getElementName());
			newElement.addText(dataFiles[i - 1].getName());
			dataFilesIDMap.put(dataFiles[i - 1], i);
		}

		// <ROW>
		int numOfRows = peakList.getNumberOfRows();
		PeakListRow row;
		for (int i = 0; i < numOfRows; i++) {
			row = peakList.getRow(i);
			newElement = saveRoot.addElement(PeakListElementName.ROW.getElementName());
			fillRowElement(row, newElement);
		}

		zipOutputStream.putNextEntry(new ZipEntry(peakList.getName()));
		OutputStream finalStream = zipOutputStream;
		OutputFormat format = OutputFormat.createPrettyPrint();
		XMLWriter writer = new XMLWriter(finalStream, format);
		writer.write(document);
	}

	
	/**
	 * @param row
	 * @param element
	 */
	private void fillRowElement(PeakListRow row, Element element) {
		element.addAttribute(PeakListElementName.ID.getElementName(), String.valueOf(row.getID()));
		Element newElement;

		// <PEAK_IDENTITY>
		PeakIdentity preferredIdentity = row.getPreferredPeakIdentity();
		PeakIdentity[] identities = row.getPeakIdentities();

		for (int i = 0; i < identities.length; i++) {
			newElement = element.addElement(PeakListElementName.PEAK_IDENTITY.getElementName());
			newElement.addAttribute(PeakListElementName.ID.getElementName(),
					String.valueOf(i));
			newElement.addAttribute(PeakListElementName.PREFERRED.getElementName(), String.valueOf(identities[i] == preferredIdentity));
			fillIdentityElement(identities[i], newElement);
		}

		// <PEAK>
		int dataFileID = 0;
		ChromatographicPeak[] peaks = row.getPeaks();
		for (ChromatographicPeak p : peaks) {
			newElement = element.addElement(PeakListElementName.PEAK.getElementName());
			dataFileID = dataFilesIDMap.get(p.getDataFile());
			fillPeakElement(p, newElement, dataFileID);
		}

	}

	/**
	 * @param identity
	 * @param element
	 */
	private void fillIdentityElement(PeakIdentity identity, Element element) {

		// <NAME>
		Element newElement = element.addElement(PeakListElementName.NAME.getElementName());
		newElement.addText(identity.getName() != null ? identity.getName() : " ");

		// <FORMULA>
		newElement = element.addElement(PeakListElementName.FORMULA.getElementName());
		String formula = "";
		if (identity instanceof SimplePeakIdentity) {
			SimplePeakIdentity id = (SimplePeakIdentity) identity;
			formula = id.getCompoundFormula();
		}
		newElement.addText(formula);

		// <IDENTIFICATION>
		newElement = element.addElement(PeakListElementName.IDENTIFICATION.getElementName());
		newElement.addText(identity.getIdentificationMethod() != null ? identity.getIdentificationMethod() : " ");

	}

	/**
	 * @param peak
	 * @param element
	 * @param dataFileID
	 */
	private void fillPeakElement(ChromatographicPeak peak, Element element,
			int dataFileID) {

		element.addAttribute(PeakListElementName.COLUMN.getElementName(),
				String.valueOf(dataFileID));
		element.addAttribute(PeakListElementName.MASS.getElementName(), String.valueOf(peak.getMZ()));
		element.addAttribute(PeakListElementName.RT.getElementName(), String.valueOf(peak.getRT()));
		element.addAttribute(PeakListElementName.HEIGHT.getElementName(),
				String.valueOf(peak.getHeight()));
		element.addAttribute(PeakListElementName.AREA.getElementName(), String.valueOf(peak.getArea()));
		element.addAttribute(PeakListElementName.STATUS.getElementName(), peak.getPeakStatus().toString());

		// <MZPEAK>
		int scanNumbers[] = peak.getScanNumbers();
		Element newElement = element.addElement(PeakListElementName.MZPEAK.getElementName());
		newElement.addAttribute(PeakListElementName.QUANTITY.getElementName(),
				String.valueOf(scanNumbers.length));


		ByteArrayOutputStream byteScanStream = new ByteArrayOutputStream();
		DataOutputStream dataScanStream = new DataOutputStream(byteScanStream);


		for (int scan : scanNumbers) {
			try {
				dataScanStream.writeInt(scan);
				dataScanStream.flush();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		char[] bytes = Base64.encode(byteScanStream.toByteArray());
		Element secondNewElement = newElement.addElement(PeakListElementName.SCAN_ID.getElementName());
		secondNewElement.addText(new String(bytes));
	}



















	public void readPeakList() {
		SAXParserFactory factory = SAXParserFactory.newInstance();
		try {
			zipInputStream.getNextEntry();
			buildingArrayRawDataFiles = new TreeMap<Integer, RawDataFile>();
			charBuffer = new StringBuffer();
			appliedProcess = new Vector<String>();
			SAXParser saxParser = factory.newSAXParser();
			saxParser.parse(zipInputStream, this);

		} catch (Throwable e) {
			e.printStackTrace();
			return;
		}

		if (parsedRows == 0) {
			return;
		}

		// Add new peaklist to the project or MZviewer.desktop
		MZmineProject currentProject = MZmineCore.getCurrentProject();
		currentProject.addPeakList(buildingPeakList);


	}

	/**
	 * @see org.xml.sax.helpers.DefaultHandler#startElement(java.lang.String,
	 *      java.lang.String, java.lang.String, org.xml.sax.Attributes)
	 */
	public void startElement(String namespaceURI, String lName, // local name
			String qName, // qualified name
			Attributes attrs) throws SAXException {

		// <PEAKLIST>
		if (qName.equals(PeakListElementName.PEAKLIST.getElementName())) {
			peakListFlag = true;
		// clean the current char buffer for the new element
		}

		// <RAWFILE>
		if (qName.equals(PeakListElementName.RAWFILE.getElementName())) {
			try {
				rawDataFileID = Integer.parseInt(attrs.getValue(PeakListElementName.ID.getElementName()));
				peakListFlag = false;
			} catch (Exception e) {
				throw new SAXException(
						"Could not read scan attributes information");
			}
		}

		// <SCAN>
		if (qName.equals(PeakListElementName.SCAN.getElementName())) {
			try {
				quantity = Integer.parseInt(attrs.getValue(PeakListElementName.QUANTITY.getElementName()));

				scanFlag = true;

			} catch (Exception e) {
				throw new SAXException(
						"Could not read scan attributes information");
			}
		}

		// <ROW>
		if (qName.equals(PeakListElementName.ROW.getElementName())) {

			if (buildingPeakList == null) {
				initializePeakList();
			}

			try {
				int rowID = Integer.parseInt(attrs.getValue(PeakListElementName.ID.getElementName()));
				buildingRow = new SimplePeakListRow(rowID);
			} catch (Exception e) {
				throw new SAXException(
						"Could not read row attributes information");
			}
		}

		// <PEAK_IDENTITY>
		if (qName.equals(PeakListElementName.PEAK_IDENTITY.getElementName())) {
			try {
				identityID = attrs.getValue(PeakListElementName.ID.getElementName());
				preferred = Boolean.parseBoolean(attrs.getValue(PeakListElementName.PREFERRED.getElementName()));
			} catch (Exception e) {
				throw new SAXException(
						"Could not read identity attributes information");
			}
		}

		// <PEAK>
		if (qName.equals(PeakListElementName.PEAK.getElementName())) {
			try {

				peakColumnID = Integer.parseInt(attrs.getValue(PeakListElementName.COLUMN.getElementName()));
				mass = Double.parseDouble(attrs.getValue(PeakListElementName.MASS.getElementName()));
				rt = Double.parseDouble(attrs.getValue(PeakListElementName.RT.getElementName()));
				height = Double.parseDouble(attrs.getValue(PeakListElementName.HEIGHT.getElementName()));
				area = Double.parseDouble(attrs.getValue(PeakListElementName.AREA.getElementName()));
				peakStatus = attrs.getValue(PeakListElementName.STATUS.getElementName());
			} catch (Exception e) {
				throw new SAXException(
						"Could not read peak attributes information");
			}

		}

		// <MZPEAK>
		if (qName.equals(PeakListElementName.MZPEAK.getElementName())) {
			try {
				quantity = Integer.parseInt(attrs.getValue(PeakListElementName.QUANTITY.getElementName()));

				mzPeakFlag = true;

			} catch (Exception e) {
				throw new SAXException(
						"Could not read mzPeak attributes information");
			}

		}

	}

	/**
	 * @see org.xml.sax.helpers.DefaultHandler#endElement(java.lang.String,
	 *      java.lang.String, java.lang.String)
	 */
	public void endElement(String namespaceURI, String sName, // simple name
			String qName // qualified name
			) throws SAXException {


		// <NAME>
		if (qName.equals(PeakListElementName.NAME.getElementName())) {
			name = getTextOfElement();
			if (peakListFlag) {
				peakListName = name;
			}

		}

		// <PEAKLIST_DATE>
		if (qName.equals(PeakListElementName.PEAKLIST_DATE.getElementName())) {
			try {
				// String text = getTextOfElement();
				dateCreated = getTextOfElement();
			} catch (Exception e) {
				throw new SAXException(
						"Could not read peak list date of creation");
			}
		}

		// <QUANTITY>
		if (qName.equals(PeakListElementName.QUANTITY.getElementName())) {
			try {
				String text = getTextOfElement();
				text = text.trim();
			//totalRows = Integer.parseInt(text);
			} catch (Exception e) {
				throw new SAXException("Could not read quantity");
			}
		}

		// <PROCESS>
		if (qName.equals(PeakListElementName.PROCESS.getElementName())) {
			String text = getTextOfElement();
			if (text.length() != 0) {
				appliedProcess.add(text);
			}
		}

		// <SCAN_ID>
		if (qName.equals(PeakListElementName.SCAN_ID.getElementName())) {
			try {
				if (scanFlag) {
					String valueText = getTextOfElement();
					String values[] = valueText.split(PeakListElementName.SEPARATOR.getElementName());
					scanNumbers = new int[quantity];
					for (int i = 0; i < quantity; i++) {
						scanNumbers[i] = Integer.parseInt(values[i]);
					}
				} else if (mzPeakFlag) {
					byte[] bytes = Base64.decode(getTextOfElement().toCharArray());
					// make a data input stream
					DataInputStream dataInputStream = new DataInputStream(
							new ByteArrayInputStream(bytes));
					scanNumbers = new int[quantity];
					for (int i = 0; i < quantity; i++) {
						scanNumbers[i] = dataInputStream.readInt();
					}
				}

			} catch (Exception e) {
				e.printStackTrace();
				throw new SAXException("Could not read list of scan numbers");
			}
		}

		// <RT>
		if (qName.equals(PeakListElementName.RT.getElementName())) {
			try {
				String valueText = getTextOfElement();
				String values[] = valueText.split(PeakListElementName.SEPARATOR.getElementName());
				retentionTimes = new double[quantity];
				for (int i = 0; i < quantity; i++) {
					retentionTimes[i] = Double.parseDouble(values[i]);
				}
			} catch (Exception e) {
				e.printStackTrace();
				throw new SAXException("Could not read list of retention times");
			}
		}

		// <MASS>
		if (qName.equals(PeakListElementName.MASS.getElementName())) {
			try {
				byte[] bytes = Base64.decode(getTextOfElement().toCharArray());
				// make a data input stream
				DataInputStream dataInputStream = new DataInputStream(
						new ByteArrayInputStream(bytes));
				masses = new double[quantity];
				for (int i = 0; i < quantity; i++) {
					masses[i] = (double) dataInputStream.readFloat();
				}

			} catch (Exception e) {
				e.printStackTrace();
				throw new SAXException("Could not read list of masses");
			}

		}

		// <HEIGHT>
		if (qName.equals(PeakListElementName.HEIGHT.getElementName())) {
			try {
				byte[] bytes = Base64.decode(getTextOfElement().toCharArray());
				// make a data input stream
				DataInputStream dataInputStream = new DataInputStream(
						new ByteArrayInputStream(bytes));
				intensities = new double[quantity];
				for (int i = 0; i < quantity; i++) {
					intensities[i] = (double) dataInputStream.readFloat();
				}

			} catch (Exception e) {
				e.printStackTrace();
				throw new SAXException("Could not read list of intensities");
			}

		}

		// <FORMULA>
		if (qName.equals(PeakListElementName.FORMULA.getElementName())) {
			formula = getTextOfElement();
		}

		// <IDENTIFICATION>
		if (qName.equals(PeakListElementName.IDENTIFICATION.getElementName())) {
			identificationMethod = getTextOfElement();
		}

		// <RTRANGE>
		if (qName.equals(PeakListElementName.RTRANGE.getElementName())) {
			try {
				String valueText = getTextOfElement();
				String values[] = valueText.split("-");
				double min = Double.parseDouble(values[0]);
				double max = Double.parseDouble(values[1]);
				rtRange = new Range(min, max);
			} catch (Exception e) {
				throw new SAXException(
						"Could not read retention time range form raw data file");
			}
		}

		// <MZRANGE>
		if (qName.equals(PeakListElementName.MZRANGE.getElementName())) {
			try {
				String valueText = getTextOfElement();
				String values[] = valueText.split("-");
				double min = Double.parseDouble(values[0]);
				double max = Double.parseDouble(values[1]);
				mzRange = new Range(min, max);
			} catch (Exception e) {
				throw new SAXException(
						"Could not read m/z range from raw data file");
			}
		}

		// <MZPEAK>
		if (qName.equals(PeakListElementName.MZPEAK.getElementName())) {
			mzPeakFlag = false;
		}

		// <PEAK>
		if (qName.equals(PeakListElementName.PEAK.getElementName())) {

			DataPoint[] mzPeaks = new DataPoint[quantity];
			Range peakRTRange = null, peakMZRange = null, peakIntensityRange = null;
			for (int i = 0; i < quantity; i++) {
				double rt = buildingArrayRawDataFiles.get(peakColumnID).getScan(scanNumbers[i]).getRetentionTime();
				double mz = masses[i];
				double intensity = intensities[i];
				if (i == 0) {
					peakRTRange = new Range(rt);
					peakMZRange = new Range(mz);
					peakIntensityRange = new Range(intensity);
				} else {
					peakRTRange.extendRange(rt);
					peakMZRange.extendRange(mz);
					peakIntensityRange.extendRange(intensity);
				}

				mzPeaks[i] = new SimpleDataPoint(mz, intensity);
			}

			SimpleChromatographicPeak peak = new SimpleChromatographicPeak(
					buildingArrayRawDataFiles.get(peakColumnID), mass, rt,
					height, area, scanNumbers, mzPeaks, PeakStatus.valueOf(
					PeakStatus.class, peakStatus), -1, -1, peakRTRange,
					peakMZRange, peakIntensityRange);

			buildingRow.addPeak(buildingArrayRawDataFiles.get(peakColumnID),
					peak);
		}

		// <PEAK_IDENTITY>
		if (qName.equals(PeakListElementName.PEAK_IDENTITY.getElementName())) {
			SimplePeakIdentity identity = new SimplePeakIdentity(identityID,
					name, new String[0], formula, null, identificationMethod);
			buildingRow.addPeakIdentity(identity, preferred);
		}

		// <ROW>
		if (qName.equals(PeakListElementName.ROW.getElementName())) {
			buildingPeakList.addRow(buildingRow);
			buildingRow = null;
			parsedRows++;
		}

		// <SCAN>
		if (qName.equals(PeakListElementName.SCAN.getElementName())) {
			try {
				if (buildingRawDataFile == null) {
					buildingRawDataFile = new RawDataFileImpl(name);
				}

				for (int i = 0; i < quantity; i++) {
					SimpleScan newScan = new SimpleScan(buildingRawDataFile,
							scanNumbers[i], 1, retentionTimes[i], -1, 0f, null,
							new DataPoint[0], false);
					buildingRawDataFile.addScan(newScan);
				}
				scanFlag = false;

			} catch (Exception e) {
				throw new SAXException(
						"Could not create scans for temporary raw data file");
			}
		}

		RawDataFile NewRawDataFile = buildingRawDataFile;

		for (RawDataFile rawDataFile : MZmineCore.getCurrentProject().getDataFiles()) {

			if (buildingRawDataFile != null && rawDataFile.getName().matches(buildingRawDataFile.getName())) {
				NewRawDataFile = (RawDataFileImpl) rawDataFile;
			}
			break;

		}

		// <RAWFILE>
		if (qName.equals(PeakListElementName.RAWFILE.getElementName())) {
			buildingRawDataFile.setRTRange(1, rtRange);
			buildingRawDataFile.setMZRange(1, mzRange);
			buildingArrayRawDataFiles.put(rawDataFileID, NewRawDataFile);
			buildingRawDataFile = null;
		}

	}

	/**
	 * Return a string without tab an EOF characters
	 *
	 * @return String element text
	 */
	private String getTextOfElement() {
		String text = charBuffer.toString();
		text = text.replaceAll("[\n\r\t]+", "");
		text = text.replaceAll("^\\s+", "");
		charBuffer.delete(0, charBuffer.length());
		return text;
	}

	/**
	 * characters()
	 *
	 * @see org.xml.sax.ContentHandler#characters(char[], int, int)
	 */
	public void characters(char buf[], int offset, int len) throws SAXException {
		charBuffer = charBuffer.append(buf, offset, len);
	}

	private void initializePeakList() {
		RawDataFile[] dataFiles = buildingArrayRawDataFiles.values().toArray(
				new RawDataFile[0]);
		buildingPeakList = new SimplePeakList(peakListName, dataFiles);
		String[] process = appliedProcess.toArray(new String[0]);
		for (String description : process) {
			((SimplePeakList) buildingPeakList).addDescriptionOfAppliedTask(new SimplePeakListAppliedMethod(
					description));
		}
		((SimplePeakList) buildingPeakList).setDateCreated(dateCreated);
	}
}