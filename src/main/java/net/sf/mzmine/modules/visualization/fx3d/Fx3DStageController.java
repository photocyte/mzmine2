/*
 * Copyright 2006-2018 The MZmine 2 Development Team
 * 
 * This file is part of MZmine 2.
 * 
 * MZmine 2 is free software; you can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 * 
 * MZmine 2 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
 * even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with MZmine 2; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
 * USA
 */
package net.sf.mzmine.modules.visualization.fx3d;

import java.util.ArrayList;

import com.google.common.collect.Range;

import javafx.fxml.FXML;
import javafx.scene.Group;
import javafx.scene.input.MouseEvent;
import javafx.scene.input.ScrollEvent;
import javafx.scene.layout.StackPane;
import javafx.scene.paint.Color;
import javafx.scene.transform.Rotate;
import javafx.scene.transform.Translate;

public class Fx3DStageController {

    @FXML
    private StackPane root = new StackPane();
    @FXML
    private Group plot = new Group();
    @FXML
    private Group finalNode = new Group();
    @FXML
    private Fx3DAxes axes = new Fx3DAxes();

    private final Rotate rotateX = new Rotate(45, Rotate.X_AXIS);
    private final Rotate rotateY = new Rotate(0, Rotate.Y_AXIS);
    private final Translate translateX = new Translate();
    private final Translate translateY = new Translate();

    private double mousePosX, mousePosY;
    private double mouseOldX, mouseOldY;

    final double MAX_SCALE = 20.0;
    final double MIN_SCALE = 0.1;

    public double maxOfAllBinnedIntensity = Double.NEGATIVE_INFINITY;

    public ArrayList<Fx3DDataset> datasets = new ArrayList<Fx3DDataset>();
    public ArrayList<Color> colors = new ArrayList<Color>();

    public void initialize() {
        plot.getTransforms().addAll(rotateX, rotateY);
        finalNode.getTransforms().addAll(translateX, translateY);
        colors.add(Color.BLUE);
        colors.add(Color.GREEN);
        colors.add(Color.RED);
        colors.add(Color.YELLOW);
        colors.add(Color.DARKORANGE);
        colors.add(Color.CYAN);
        colors.add(Color.FUCHSIA);
        colors.add(Color.GOLD);
    }

    public synchronized void setDataset(Fx3DDataset dataset,
            double maxBinnedIntensity, int index, int length) {
        datasets.add(dataset);
        if (maxOfAllBinnedIntensity < maxBinnedIntensity) {
            maxOfAllBinnedIntensity = maxBinnedIntensity;
        }
        if (index == length - 1) {
            int i = 0;
            for (Fx3DDataset data : datasets) {
                Fx3DPlotMesh meshView = new Fx3DPlotMesh();
                meshView.setDataset(data, maxOfAllBinnedIntensity,
                        colors.get(i));
                plot.getChildren().addAll(meshView);
                i = (i + 1) % 8;
            }
            Range<Double> rtRange = dataset.getRtRange();
            Range<Double> mzRange = dataset.getMzRange();
            axes.setValues(rtRange, mzRange, maxOfAllBinnedIntensity);
        }
    }

    public void handleMousePressed(MouseEvent me) {
        mouseOldX = me.getSceneX();
        mouseOldY = me.getSceneY();
    }

    public void handleMouseDragged(MouseEvent me) {
        double rotateFactor = 0.08;
        mousePosX = me.getSceneX();
        mousePosY = me.getSceneY();
        if (me.isPrimaryButtonDown()) {
            rotateX.setAngle(rotateX.getAngle()
                    + rotateFactor * (mousePosY - mouseOldY));
            rotateY.setAngle(rotateY.getAngle()
                    - rotateFactor * (mousePosX - mouseOldX));
        }
        if (me.isSecondaryButtonDown()) {
            translateX.setX(translateX.getX() + (mousePosX - mouseOldX));
            translateY.setY(translateY.getY() + (mousePosY - mouseOldY));
        }
        mouseOldX = mousePosX;
        mouseOldY = mousePosY;
    }

    public void onScrollHandler(ScrollEvent event) {
        double delta = 1.2;
        double scale = (root.getScaleX());

        if (event.getDeltaY() < 0) {
            scale /= delta;
        } else {
            scale *= delta;
        }

        scale = clamp(scale, MIN_SCALE, MAX_SCALE);

        root.setScaleX(scale);
        root.setScaleY(scale);

        event.consume();
    }

    public static double clamp(double value, double min, double max) {

        if (Double.compare(value, min) < 0)
            return min;

        if (Double.compare(value, max) > 0)
            return max;

        return value;
    }

    public Fx3DAxes getAxes() {
        return axes;
    }
}
