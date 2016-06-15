import os
import logging

from astropy.units import Unit

from ..core.comms import Dispatch, DispatchHandle
from ..ui.widgets.utils import ICON_PATH
from ..ui.widgets.plugin import Plugin
from ..ui.widgets.dialogs import TopAxisDialog, UnitChangeDialog, LineListsWindow


class PlotToolsPlugin(Plugin):
    name = "Plot Tools"
    location = "hidden"
    _all_categories = {}

    def setup_ui(self):
        self._top_axis_dialog = TopAxisDialog()
        self._unit_change_dialog = UnitChangeDialog()
        self._linelist_window = LineListsWindow()

        # Add an roi
        self.button_add_roi = self.add_tool_bar_actions(
            name="ROI",
            description='Add ROI',
            icon_path=os.path.join(ICON_PATH, "Merge Vertical-48.png"),
            category='Selections',
            callback=Dispatch.on_add_roi.emit,
            enabled=False)

        # Change top axis
        self.button_axis_change = self.add_tool_bar_actions(
            name="Top Axis",
            description='Change top axis',
            icon_path=os.path.join(ICON_PATH, "Globe Earth-48.png"),
            category='Options',
            callback=self._top_axis_dialog.exec_,
            enabled=False)

        # Change top axis
        self.button_unit_change = self.add_tool_bar_actions(
            name="Units",
            description='Change plot units',
            icon_path=os.path.join(ICON_PATH, "Generic Text-48.png"),
            category='Options',
            callback=self._show_unit_change_dialog,
            enabled=False)

        self.button_line_labels = self.add_tool_bar_actions(
            name="Line Labels",
            description='Add line labels',
            icon_path=os.path.join(ICON_PATH, "Label-48.png"),
            category='Selections',
            callback=self._show_linelists_window,
            enabled=False)

    def setup_connections(self):
        # On accept, change the displayed axis
        self._top_axis_dialog.accepted.connect(
            self._update_axis)

    def _show_unit_change_dialog(self):
        if self._unit_change_dialog.exec_():
            x_text = self._unit_change_dialog.disp_unit
            y_text = self._unit_change_dialog.flux_unit

            x_unit = y_unit = None

            try:
                x_unit = Unit(x_text) if x_text else None
            except ValueError as e:
                logging.error(e)

            try:
                y_unit = Unit(y_text) if y_text else None
            except ValueError as e:
                logging.error(e)

            self.active_window.change_units(x_unit, y_unit)
            self.active_window.update_plot_item()

    def _update_axis(self):
        if self.active_window is None:
            return

        if len(self.active_window._plots) > 0:
            layer = self.active_window._plots[0].layer

            self.active_window.update_axis(
                layer,
                self._top_axis_dialog.combo_box_axis_mode.currentIndex(),
                redshift=self._top_axis_dialog.redshift,
                ref_wave=self._top_axis_dialog.ref_wave)
        else:
            logging.warning("Active window does not have any plots.")

#TODO work in progress. Behavior when closing, TBD.
    def _show_linelists_window(self):
        self._linelist_window.show()

    @DispatchHandle.register_listener("on_activated_window")
    def toggle_enabled(self, window):
        if window:
            self.button_axis_change.setEnabled(True)
            self.button_unit_change.setEnabled(True)
            self.button_add_roi.setEnabled(True)
            self.button_line_labels.setEnabled(True)
        else:
            self.button_axis_change.setEnabled(False)
            self.button_unit_change.setEnabled(False)
            self.button_add_roi.setEnabled(False)
            self.button_line_labels.setEnabled(False)
