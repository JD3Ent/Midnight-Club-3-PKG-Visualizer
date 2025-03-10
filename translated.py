import sys
import struct
import random
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QFileDialog, QMessageBox, QAction,
    QDockWidget, QTreeWidget, QTreeWidgetItem, QWidget, QVBoxLayout, QCheckBox,
    QPushButton, QColorDialog, QLabel, QHBoxLayout, QSlider, QComboBox
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QColor
import pyvista as pv
from pyvistaqt import QtInteractor
import numpy as np


class PCKParser:
    """
    Class to read and parse PCK files according to provided specification.
    Now UVs are read and associated with each group of vertices.
    """

    def __init__(self, filepath):
        self.filepath = filepath
        self.vertex_groups = []  # List of groups with vertices, faces and UVs
        self.parse_success = False
        self.error_message = ""
        self.data = b''

    def read_pck_file(self):
        try:
            with open(self.filepath, 'rb') as f:
                self.data = f.read()
            return True
        except Exception as e:
            self.error_message = f"Error reading the file: {e}"
            return False

    def find_vertex_patterns(self):
        """
        Finds vertex and UV patterns in the file.

        Vertex Patterns:
          - pattern1: EE 00 XX 69
          - pattern2: 1B 02 XX 69
        Where XX is between 00 and 2A and each vertex occupies 6 bytes.

        After vertices, we try to read UVs:
          - For pattern1, expected header: C4 00 XX 65
          - For pattern2, expected header: F1 01 XX 65
        Each UV pair has 4 bytes (2 bytes for U and 2 bytes for V, 16-bit signed).

        Next, we look for the face pattern:
          - For pattern1: 9A00 + XX + 6A
          - For pattern2: C701 + XX + 6A
        """
        patterns = {
            'pattern1': {
                'prefix': bytes.fromhex('EE00'),
                'faces_prefix': bytes.fromhex('9A00'),
                'uv_prefix': bytes.fromhex('C400'),
                'type': 'vertices'
            },
            'pattern2': {
                'prefix': bytes.fromhex('1B02'),
                'faces_prefix': bytes.fromhex('C701'),
                'uv_prefix': bytes.fromhex('F101'),
                'type': 'vertices'
            }
        }
        i = 0
        while i < len(self.data) - 4:
            found_pattern = False
            for pattern_name, pattern_info in patterns.items():
                prefix = pattern_info['prefix']
                if self.data[i:i+2] == prefix:
                    XX = self.data[i+2]
                    if not (0x00 <= XX <= 0x2A):
                        continue  # XX outside expected range
                    # Check the final byte of the vertex header (must be 0x69)
                    if self.data[i+3] != 0x69:
                        continue

                    found_pattern = True

                    group = {
                        'pattern': pattern_name,
                        'start': i,
                        'XX': XX,
                        'vertices': [],
                        'active_faces': [],
                        'inactive_faces': [],
                        'uvs': [],  # List of pairs (U,V)
                        'color': self.generate_group_color()
                    }
                    # Read vertices
                    vertices_start = i + 4
                    vertex_bytes_count = XX * 6  # Each vertex occupies 6 bytes
                    if vertices_start + vertex_bytes_count > len(self.data):
                        self.error_message = "Corrupted or incomplete PCK file (vertices)."
                        return False
                    for j in range(0, vertex_bytes_count, 6):
                        try:
                            x, y, z = struct.unpack_from('<hhh', self.data, vertices_start + j)
                            group['vertices'].append((x, y, z))
                        except struct.error as e:
                            self.error_message = f"Error parsing vertices: {e}"
                            return False
                    group['hex_vertices'] = self.data[vertices_start:vertices_start + vertex_bytes_count].hex().upper()
                    group['vertices_offset'] = vertices_start
                    group['vertices_end'] = vertices_start + vertex_bytes_count

                    # Read UVs immediately after vertices, if there's an expected header
                    pos_after_vertices = group['vertices_end']
                    if pos_after_vertices + 4 <= len(self.data):
                        expected_uv_header = pattern_info['uv_prefix'] + bytes([XX]) + bytes.fromhex('65')
                        if self.data[pos_after_vertices: pos_after_vertices+4] == expected_uv_header:
                            uv_start = pos_after_vertices + 4
                            uv_bytes_count = XX * 4  # Each UV pair occupies 4 bytes
                            if uv_start + uv_bytes_count > len(self.data):
                                self.error_message = "Corrupted or incomplete PCK file (UVs)."
                                return False
                            for j in range(0, uv_bytes_count, 4):
                                try:
                                    u, v = struct.unpack_from('<hh', self.data, uv_start + j)
                                    group['uvs'].append((u, v))
                                except struct.error as e:
                                    self.error_message = f"Error parsing UVs: {e}"
                                    return False
                            group['hex_uvs'] = self.data[uv_start: uv_start + uv_bytes_count].hex().upper()
                            group['uvs_offset'] = uv_start
                            group['uvs_end'] = uv_start + uv_bytes_count
                            pos_after_uv = group['uvs_end']
                        else:
                            group['uvs'] = []
                            pos_after_uv = pos_after_vertices
                    else:
                        group['uvs'] = []
                        pos_after_uv = pos_after_vertices

                    # Search for face pattern starting from pos_after_uv
                    faces_prefix = pattern_info['faces_prefix'] + bytes([XX]) + bytes.fromhex('6A')
                    pos_faces = self.data.find(faces_prefix, pos_after_uv)
                    if pos_faces == -1:
                        self.error_message = "Face pattern not found."
                        return False
                    faces_start = pos_faces + len(faces_prefix) + 6  # Skip 6 bytes after face header
                    faces_bytes_count = (XX - 2) * 3
                    if faces_start + faces_bytes_count > len(self.data):
                        self.error_message = "Corrupted or incomplete PCK file (faces)."
                        return False
                    for j in range(0, faces_bytes_count, 3):
                        try:
                            face_bytes = struct.unpack_from('<BBB', self.data, faces_start + j)
                            value = face_bytes[0]
                            face_idx = j // 3 + 1
                            v1 = face_idx
                            v2 = face_idx + 1
                            v3 = face_idx + 2
                            if v3 > XX:
                                continue
                            face = {'indices': (v1, v2, v3), 'byte': value}
                            if value % 2 == 0:
                                group['active_faces'].append(face)
                            else:
                                group['inactive_faces'].append(face)
                        except struct.error as e:
                            self.error_message = f"Error parsing faces: {e}"
                            return False
                    group['hex_faces'] = self.data[pos_faces: pos_faces + len(faces_prefix) + faces_bytes_count].hex().upper()
                    group['faces_offset'] = pos_faces
                    group['faces_end'] = pos_faces + len(faces_prefix) + faces_bytes_count

                    self.vertex_groups.append(group)
                    i = faces_start + faces_bytes_count
                    break
            if not found_pattern:
                i += 1

        print(f"Total Vertex Groups Found: {len(self.vertex_groups)}")
        for idx, group in enumerate(self.vertex_groups):
            print(f"Group {idx + 1}:")
            print(f"  Pattern: {group['pattern']}")
            print(f"  Vertex Count: {group['XX']}")
            print(f"  Vertex Offset: {hex(group['vertices_offset'])} - {hex(group['vertices_end'])}")
            print(f"  UVs Found: {len(group['uvs'])}")
            print(f"  Active Faces: {len(group['active_faces'])}")
            print(f"  Inactive Faces: {len(group['inactive_faces'])}")

        total_vertices = sum(len(group['vertices']) for group in self.vertex_groups)
        total_active_faces = sum(len(group['active_faces']) for group in self.vertex_groups)
        total_inactive_faces = sum(len(group['inactive_faces']) for group in self.vertex_groups)
        print(f"Total Global Vertices: {total_vertices}")
        print(f"Total Global Active Faces: {total_active_faces}")
        print(f"Total Global Inactive Faces: {total_inactive_faces}")

        self.parse_success = True
        return True

    def generate_group_color(self):
        return (random.random(), random.random(), random.random(), 1)

    def parse(self):
        if not self.read_pck_file():
            return False
        success = self.find_vertex_patterns()
        if not success:
            return False

        # Adjust global vertex indices for faces
        self.active_faces = []
        self.inactive_faces = []
        self.vertices = []
        for group in self.vertex_groups:
            base_index = len(self.vertices)
            self.vertices.extend(group['vertices'])
            group['base_index'] = base_index
            group['mesh_active_faces'] = []
            group['mesh_inactive_faces'] = []
            for face in group['active_faces']:
                if isinstance(face, dict) and 'indices' in face:
                    adjusted_indices = tuple(idx - 1 + base_index for idx in face['indices'])
                    if all(0 <= idx < len(self.vertices) for idx in adjusted_indices):
                        group['mesh_active_faces'].append(adjusted_indices)
                else:
                    self.error_message = f"Invalid face format in group {group['pattern']}."
                    return False
            for face in group['inactive_faces']:
                if isinstance(face, dict) and 'indices' in face:
                    adjusted_indices = tuple(idx - 1 + base_index for idx in face['indices'])
                    if all(0 <= idx < len(self.vertices) for idx in adjusted_indices):
                        group['mesh_inactive_faces'].append(adjusted_indices)
                else:
                    self.error_message = f"Invalid face format in group {group['pattern']}."
                    return False

        print(f"Total Global Vertices: {len(self.vertices)}")
        total_active_faces = sum(len(group['mesh_active_faces']) for group in self.vertex_groups)
        total_inactive_faces = sum(len(group['mesh_inactive_faces']) for group in self.vertex_groups)
        print(f"Total Global Active Faces: {total_active_faces}")
        print(f"Total Global Inactive Faces: {total_inactive_faces}")

        self.parse_success = True
        return True


class InformationPanel(QDockWidget):
    def __init__(self, parent=None):
        super(InformationPanel, self).__init__("Information", parent)
        self.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(['Property', 'Value'])
        self.setWidget(self.tree)

    def update_info(self, vertex_groups):
        self.tree.clear()
        for idx, group in enumerate(vertex_groups):
            group_item = QTreeWidgetItem(self.tree, [f'Group {idx + 1} - {group["pattern"]}', ''])
            color = group['color']
            qcolor = QColor(int(color[0]*255), int(color[1]*255), int(color[2]*255), int(color[3]*255))
            group_item.setBackground(0, qcolor)
            vertices_item = QTreeWidgetItem(group_item, ['Vertices', ''])
            QTreeWidgetItem(vertices_item, ['Count', str(group['XX'])])
            QTreeWidgetItem(vertices_item, ['Size (bytes)', str(group['vertices_end'] - group['vertices_offset'])])
            QTreeWidgetItem(vertices_item, ['Start Offset', hex(group['vertices_offset'])])
            QTreeWidgetItem(vertices_item, ['End Offset', hex(group['vertices_end'])])
            QTreeWidgetItem(vertices_item, ['Hex Vertices', group.get('hex_vertices', '')])
            uvs_item = QTreeWidgetItem(group_item, ['UVs', ''])
            QTreeWidgetItem(uvs_item, ['Count', str(len(group['uvs']))])
            QTreeWidgetItem(uvs_item, ['Start Offset', hex(group.get('uvs_offset', 0))])
            QTreeWidgetItem(uvs_item, ['End Offset', hex(group.get('uvs_end', 0))])
            QTreeWidgetItem(uvs_item, ['Hex UVs', group.get('hex_uvs', '')])
            active_faces_item = QTreeWidgetItem(group_item, ['Active Faces', ''])
            QTreeWidgetItem(active_faces_item, ['Count', str(len(group['mesh_active_faces']))])
            for face_idx, face in enumerate(group['mesh_active_faces']):
                face_info = f'Face {face_idx + 1}: {face}'
                QTreeWidgetItem(active_faces_item, [face_info, 'Active'])
            inactive_faces_item = QTreeWidgetItem(group_item, ['Inactive Faces', ''])
            QTreeWidgetItem(inactive_faces_item, ['Count', str(len(group['mesh_inactive_faces']))])
            for face_idx, face in enumerate(group['mesh_inactive_faces']):
                face_info = f'Face {face_idx + 1}: {face}'
                QTreeWidgetItem(inactive_faces_item, [face_info, 'Inactive'])
            total_faces = len(group['mesh_active_faces']) + len(group['mesh_inactive_faces'])
            QTreeWidgetItem(group_item, ['Total Faces', str(total_faces)])
        self.tree.expandAll()


class TogglePanel(QDockWidget):
    def __init__(self, parent=None):
        super(TogglePanel, self).__init__("Inactive Faces Control", parent)
        self.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
        widget = QWidget()
        layout = QVBoxLayout()
        self.checkbox = QCheckBox("Show Inactive Faces")
        self.checkbox.setChecked(False)
        self.checkbox.setToolTip("Check to display inactive faces in 3D view.")
        layout.addWidget(self.checkbox)
        layout.addStretch()
        widget.setLayout(layout)
        self.setWidget(widget)


class SettingsPanel(QDockWidget):
    def __init__(self, parent=None):
        super(SettingsPanel, self).__init__("Settings", parent)
        self.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
        widget = QWidget()
        layout = QVBoxLayout()
        bg_color_layout = QHBoxLayout()
        bg_color_label = QLabel("Background Color:")
        self.bg_color_button = QPushButton("Choose Color")
        self.bg_color_button.setToolTip("Click to choose the 3D view background color.")
        bg_color_layout.addWidget(bg_color_label)
        bg_color_layout.addWidget(self.bg_color_button)
        layout.addLayout(bg_color_layout)
        self.wireframe_checkbox = QCheckBox("Show Wireframe")
        self.wireframe_checkbox.setChecked(False)
        self.wireframe_checkbox.setToolTip("Enable to view meshes as wireframe.")
        layout.addWidget(self.wireframe_checkbox)
        self.grid_checkbox = QCheckBox("Show Grid")
        self.grid_checkbox.setChecked(True)
        self.grid_checkbox.setToolTip("Enable to show the reference grid in the 3D scene.")
        layout.addWidget(self.grid_checkbox)
        self.axes_checkbox = QCheckBox("Show Axes")
        self.axes_checkbox.setChecked(True)
        self.axes_checkbox.setToolTip("Enable to show coordinate axes in the 3D scene.")
        layout.addWidget(self.axes_checkbox)
        export_layout = QHBoxLayout()
        export_label = QLabel("Export Model:")
        self.export_button = QPushButton("Export to OBJ")
        self.export_button.setToolTip("Click to export the current model to an OBJ file (including UVs), omitting inactive faces.")
        export_layout.addWidget(export_label)
        export_layout.addWidget(self.export_button)
        layout.addLayout(export_layout)
        transparency_layout = QHBoxLayout()
        transparency_label = QLabel("Mesh Transparency:")
        self.transparency_slider = QSlider(Qt.Horizontal)
        self.transparency_slider.setRange(0, 100)
        self.transparency_slider.setValue(100)
        self.transparency_slider.setToolTip("Adjust mesh transparency.")
        transparency_layout.addWidget(transparency_label)
        transparency_layout.addWidget(self.transparency_slider)
        layout.addLayout(transparency_layout)
        shading_layout = QHBoxLayout()
        shading_label = QLabel("Shading Type:")
        self.shading_combo = QComboBox()
        self.shading_combo.addItems(["Flat", "Smooth"])
        self.shading_combo.setToolTip("Select shading type for meshes.")
        shading_layout.addWidget(shading_label)
        shading_layout.addWidget(self.shading_combo)
        layout.addLayout(shading_layout)
        scale_layout = QHBoxLayout()
        scale_label = QLabel("Scale Model:")
        self.scale_checkbox = QCheckBox("Divide Vertices by 256")
        self.scale_checkbox.setToolTip("Enable to scale the model by dividing all vertex coordinates by 256.")
        scale_layout.addWidget(scale_label)
        scale_layout.addWidget(self.scale_checkbox)
        layout.addLayout(scale_layout)
        layout.addStretch()
        widget.setLayout(layout)
        self.setWidget(widget)


class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.setWindowTitle("PCK Files 3D Viewer")
        self.setGeometry(100, 100, 1200, 800)
        self.view_widget = QtInteractor(self)
        self.setCentralWidget(self.view_widget)

        # Improve 3D visualization: enable anti-aliasing and configure render window
        self.view_widget.render_window.SetMultiSamples(8)
        self.view_widget.reset_camera()

        self.grid_actor = self.add_grid()
        self.axes_actor = self.add_axes()
        self.info_panel = InformationPanel(self)
        self.addDockWidget(Qt.RightDockWidgetArea, self.info_panel)
        self.toggle_panel = TogglePanel(self)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.toggle_panel)
        self.settings_panel = SettingsPanel(self)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.settings_panel)
        self.toggle_panel.checkbox.stateChanged.connect(self.toggle_inactive_faces)
        self.settings_panel.bg_color_button.clicked.connect(self.choose_background_color)
        self.settings_panel.wireframe_checkbox.stateChanged.connect(self.toggle_wireframe)
        self.settings_panel.grid_checkbox.stateChanged.connect(self.toggle_grid)
        self.settings_panel.axes_checkbox.stateChanged.connect(self.toggle_axes)
        self.settings_panel.export_button.clicked.connect(self.export_to_obj)
        self.settings_panel.transparency_slider.valueChanged.connect(self.adjust_transparency)
        self.settings_panel.shading_combo.currentIndexChanged.connect(self.change_shading)
        self.settings_panel.scale_checkbox.stateChanged.connect(self.toggle_scale)
        self.mesh_groups = []
        self.mesh_inactive = None
        self.current_shading = "Smooth"
        self.scale_enabled = False
        self.original_vertices = []
        self.init_menu()
        self.parser = None

    def init_menu(self):
        menubar = self.menuBar()
        file_menu = menubar.addMenu('File')
        open_action = QAction('Open PCK', self)
        open_action.setShortcut('Ctrl+O')
        open_action.setStatusTip('Open PCK file')
        open_action.triggered.connect(self.open_file)
        file_menu.addAction(open_action)
        export_action = QAction('Export to OBJ', self)
        export_action.setShortcut('Ctrl+E')
        export_action.setStatusTip('Export model to OBJ')
        export_action.triggered.connect(self.export_to_obj)
        file_menu.addAction(export_action)
        exit_action = QAction('Exit', self)
        exit_action.setShortcut('Ctrl+Q')
        exit_action.setStatusTip('Exit application')
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)

    def add_grid(self):
        grid_size = 1000
        spacing = 100
        points = []
        lines_encoded = []
        for y in range(-grid_size, grid_size + spacing, spacing):
            points.append([-grid_size, y, 0])
            points.append([grid_size, y, 0])
        for x in range(-grid_size, grid_size + spacing, spacing):
            points.append([x, -grid_size, 0])
            points.append([x, grid_size, 0])
        points_np = np.array(points, dtype=np.float32)
        grid = pv.PolyData(points_np)
        n_lines = len(points_np) // 2
        for i in range(n_lines):
            lines_encoded.append(2)
            lines_encoded.append(i * 2)
            lines_encoded.append(i * 2 + 1)
        lines_np = np.array(lines_encoded, dtype=np.int32)
        grid.lines = lines_np
        grid_actor = self.view_widget.add_mesh(grid, color='gray', line_width=1, opacity=0.5, name="Grid")
        return grid_actor

    def add_axes(self):
        axes_actor = self.view_widget.add_axes()
        return axes_actor

    def open_file(self):
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly
        filepath, _ = QFileDialog.getOpenFileName(
            self,
            "Open PCK File",
            "",
            "PCK Files (*.pck);;All Files (*)",
            options=options
        )
        if filepath:
            parser = PCKParser(filepath)
            success = parser.parse()
            if success:
                self.original_vertices = parser.vertices.copy()
                self.render_model(parser.vertex_groups, parser.vertices)
                self.info_panel.update_info(parser.vertex_groups)
                self.setWindowTitle(f"PCK Files 3D Viewer - {filepath}")
                self.parser = parser
            else:
                QMessageBox.critical(self, "Error", parser.error_message)

    def render_model(self, vertex_groups, vertices):
        self.view_widget.clear()
        self.grid_actor = self.add_grid()
        self.axes_actor = self.add_axes()
        if self.scale_enabled:
            vertices_np = np.array(vertices, dtype=np.float32) / 256.0
        else:
            vertices_np = np.array(vertices, dtype=np.float32)
        self.mesh_groups = []
        for idx, group in enumerate(vertex_groups):
            if group['mesh_active_faces']:
                faces_np = np.array(group['mesh_active_faces'], dtype=np.int32)
                faces_formatted = np.hstack([np.full((faces_np.shape[0], 1), 3), faces_np]).flatten()
                mesh = pv.PolyData(vertices_np, faces_formatted)
                mesh_plot = self.view_widget.add_mesh(
                    mesh,
                    color=group['color'][:3],
                    show_edges=self.settings_panel.wireframe_checkbox.isChecked(),
                    edge_color='black',
                    name=f"Group {idx + 1}",
                    smooth_shading=(self.current_shading == "Smooth")
                )
                self.mesh_groups.append(mesh_plot)
        self.mesh_inactive = None
        inactive_faces_all = []
        for group in vertex_groups:
            if group['mesh_inactive_faces']:
                inactive_faces_all.extend(group['mesh_inactive_faces'])
        if inactive_faces_all:
            inactive_faces_np = np.array(inactive_faces_all, dtype=np.int32)
            faces_formatted = np.hstack([np.full((inactive_faces_np.shape[0], 1), 3), inactive_faces_np]).flatten()
            mesh_inactive = pv.PolyData(vertices_np, faces_formatted)
            self.mesh_inactive = self.view_widget.add_mesh(
                mesh_inactive,
                color=(1, 0, 0),
                show_edges=self.settings_panel.wireframe_checkbox.isChecked(),
                edge_color='black',
                name="Inactive Faces",
                smooth_shading=(self.current_shading == "Smooth")
            )
            self.mesh_inactive.SetVisibility(False)
        self.view_widget.reset_camera()

    def toggle_inactive_faces(self, state):
        if self.mesh_inactive:
            self.mesh_inactive.SetVisibility(state == Qt.Checked)
            self.view_widget.render()

    def choose_background_color(self):
        color = QColorDialog.getColor()
        if color.isValid():
            rgb_color = (color.red() / 255, color.green() / 255, color.blue() / 255)
            self.view_widget.set_background(rgb_color)

    def toggle_wireframe(self, state):
        show_wireframe = state == Qt.Checked
        for mesh in self.mesh_groups:
            mesh.GetProperty().SetEdgeVisibility(show_wireframe)
        if self.mesh_inactive:
            self.mesh_inactive.GetProperty().SetEdgeVisibility(show_wireframe)
        self.view_widget.render()

    def toggle_grid(self, state):
        if self.grid_actor:
            self.grid_actor.SetVisibility(state == Qt.Checked)
            self.view_widget.render()

    def toggle_axes(self, state):
        if self.axes_actor:
            self.axes_actor.SetVisibility(state == Qt.Checked)
            self.view_widget.render()

    def adjust_transparency(self, value):
        transparency = value / 100.0
        for mesh in self.mesh_groups:
            # Adjust opacity using VTK actor property
            mesh.GetProperty().SetOpacity(transparency)
        if self.mesh_inactive:
            self.mesh_inactive.GetProperty().SetOpacity(transparency)
        self.view_widget.render()

    def change_shading(self, index):
        shading_types = ["Flat", "Smooth"]
        self.current_shading = shading_types[index]
        self.render_model_with_current_settings()
        self.view_widget.render()

    def toggle_scale(self, state):
        self.scale_enabled = state == Qt.Checked
        if self.parser and self.parser.parse_success:
            self.render_model_with_current_settings()

    def render_model_with_current_settings(self):
        if not self.parser or not self.parser.parse_success:
            return
        vertices = self.original_vertices
        self.render_model(self.parser.vertex_groups, vertices)

    def export_to_obj(self):
        if not self.parser or not self.parser.parse_success:
            QMessageBox.warning(self, "Warning", "No model loaded for export.")
            return

        filepath, _ = QFileDialog.getSaveFileName(
            self,
            "Export to OBJ",
            "",
            "OBJ Files (*.obj);;All Files (*)"
        )
        if not filepath:
            return

        if self.scale_enabled:
            vertices = np.array(self.original_vertices, dtype=np.float32) / 256.0
        else:
            vertices = np.array(self.original_vertices, dtype=np.float32)

        # Build grouped UVs from each group, in vertex order
        global_uvs = []
        for group in self.parser.vertex_groups:
            if len(group['uvs']) == len(group['vertices']):
                global_uvs.extend(group['uvs'])
            else:
                global_uvs.extend([(0, 0)] * len(group['vertices']))
        if len(global_uvs) < len(vertices):
            global_uvs.extend([(0, 0)] * (len(vertices) - len(global_uvs)))
        elif len(global_uvs) > len(vertices):
            global_uvs = global_uvs[:len(vertices)]

        obj_lines = []
        for v in vertices:
            obj_lines.append(f"v {v[0]} {v[1]} {v[2]}")
        for uv in global_uvs:
            obj_lines.append(f"vt {uv[0]} {uv[1]}")
        for group in self.parser.vertex_groups:
            for face in group['mesh_active_faces']:
                v1, v2, v3 = face
                obj_lines.append(f"f {v1+1}/{v1+1} {v2+1}/{v2+1} {v3+1}/{v3+1}")
        obj_content = "\n".join(obj_lines)
        print("OBJ file content:")
        print(obj_content)
        try:
            with open(filepath, "w") as f:
                f.write(obj_content)
            QMessageBox.information(self, "Success", f"Model successfully exported to {filepath}.")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to export model: {e}")

    def change_lighting(self, index):
        pass

def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
