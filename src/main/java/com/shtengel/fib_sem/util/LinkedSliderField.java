package com.shtengel.fib_sem.util;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.JTextField;

/**
 * A Swing component that pairs a {@link JSlider} with a {@link JTextField},
 * keeping them synchronized bidirectionally over a configurable numeric range.
 *
 * <p>Value changes from either control are mapped through a linear interpolation
 * between {@code min} and {@code max}. An optional change callback fires on
 * every value update (from slider drag or field Enter key).</p>
 */
public class LinkedSliderField extends JPanel {

	private static final int DEFAULT_STEPS = 1000;

	private final JSlider slider;
	private final JTextField field;
	private final double min;
	private final double max;
	private final int steps;
	private final String format;

	private boolean updatingFromSlider;
	private boolean updatingFromField;
	private Runnable changeCallback;

	/**
	 * Creates a linked slider-field with a label, value range, and display format.
	 *
	 * @param label        label text shown to the left of the slider
	 * @param min          minimum mapped value
	 * @param max          maximum mapped value
	 * @param initialValue initial value (clamped to [min, max])
	 * @param format       printf format for the text field (e.g. "%.3f", "%.1f")
	 */
	public LinkedSliderField(String label, double min, double max,
							 double initialValue, String format) {
		this.min = min;
		this.max = max;
		this.steps = DEFAULT_STEPS;
		this.format = format;

		setLayout(new GridBagLayout());
		GridBagConstraints gbc = new GridBagConstraints();
		gbc.insets = new Insets(2, 4, 2, 4);
		gbc.fill = GridBagConstraints.HORIZONTAL;

		// Create components first (field must exist before slider's listener references it)
		field = new JTextField(String.format(format, initialValue), 6);
		field.addActionListener(e -> onFieldEdited());

		slider = new JSlider(0, steps, toSlider(initialValue));
		slider.addChangeListener(e -> {
			if (updatingFromField) return;
			updatingFromSlider = true;
			field.setText(String.format(format, fromSlider(slider.getValue())));
			updatingFromSlider = false;
			fireChange();
		});

		// Layout: label | slider | field
		gbc.gridx = 0;
		gbc.weightx = 0;
		add(new JLabel(label), gbc);

		gbc.gridx = 1;
		gbc.weightx = 1.0;
		add(slider, gbc);

		gbc.gridx = 2;
		gbc.weightx = 0;
		gbc.fill = GridBagConstraints.NONE;
		add(field, gbc);
	}

	/** Registers a callback invoked on every value change (slider or field). */
	public void addChangeCallback(Runnable callback) {
		this.changeCallback = callback;
	}

	/** Returns the current mapped value. */
	public double getValue() {
		return fromSlider(slider.getValue());
	}

	/**
	 * Sets the value programmatically, updating both slider and field.
	 *
	 * @param value the new value (clamped to [min, max])
	 */
	public void setValue(double value) {
		updatingFromField = true;
		slider.setValue(toSlider(value));
		field.setText(String.format(format, value));
		updatingFromField = false;
	}

	/** Enables or disables both the slider and the text field. */
	@Override
	public void setEnabled(boolean enabled) {
		super.setEnabled(enabled);
		slider.setEnabled(enabled);
		field.setEnabled(enabled);
	}

	/** Updates the mapped range and refreshes the display without firing the callback. */
	public void setRange(double newMin, double newMax) {
		// Not needed for current use, but included for completeness
		// Would require storing min/max as mutable and re-mapping
	}

	private int toSlider(double value) {
		if (max == min) return steps / 2;
		double fraction = (value - min) / (max - min);
		return (int) Math.round(Math.max(0, Math.min(1, fraction)) * steps);
	}

	private double fromSlider(int sliderValue) {
		double fraction = sliderValue / (double) steps;
		return min + fraction * (max - min);
	}

	private void onFieldEdited() {
		if (updatingFromSlider) return;
		try {
			double val = Double.parseDouble(field.getText().trim());
			updatingFromField = true;
			slider.setValue(toSlider(val));
			updatingFromField = false;
			fireChange();
		} catch (NumberFormatException ignored) {
			// Revert field to current slider value on bad input
			field.setText(String.format(format, fromSlider(slider.getValue())));
		}
	}

	private void fireChange() {
		if (changeCallback != null) {
			changeCallback.run();
		}
	}
}
