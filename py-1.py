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
    Classe para ler e parsear arquivos PCK conforme especificação fornecida.
    Agora os UVs são lidos e associados a cada grupo de vértices.
    """

    def __init__(self, filepath):
        self.filepath = filepath
        self.vertex_groups = []  # Lista de grupos com vértices, faces e UVs
        self.parse_success = False
        self.error_message = ""
        self.dados = b''

    def ler_arquivo_pck(self):
        try:
            with open(self.filepath, 'rb') as f:
                self.dados = f.read()
            return True
        except Exception as e:
            self.error_message = f"Erro ao ler o arquivo: {e}"
            return False

    def encontrar_padroes_vertices(self):
        """
        Encontra os padrões de vértices e UVs no arquivo.

        Padrões de Vértices:
          - padrao1: EE 00 XX 69
          - padrao2: 1B 02 XX 69
        Onde XX está entre 00 e 2A e cada vértice ocupa 6 bytes.

        Após os vértices, tenta-se ler os UVs:
          - Para padrao1, cabeçalho esperado: C4 00 XX 65
          - Para padrao2, cabeçalho esperado: F1 01 XX 65
        Cada par UV tem 4 bytes (2 bytes para U e 2 bytes para V, 16-bit signed).

        Em seguida, procura-se o padrão de faces:
          - Para padrao1: 9A00 + XX + 6A
          - Para padrao2: C701 + XX + 6A
        """
        padroes = {
            'padrao1': {
                'prefixo': bytes.fromhex('EE00'),
                'faces_prefixo': bytes.fromhex('9A00'),
                'uv_prefixo': bytes.fromhex('C400'),
                'type': 'vertices'
            },
            'padrao2': {
                'prefixo': bytes.fromhex('1B02'),
                'faces_prefixo': bytes.fromhex('C701'),
                'uv_prefixo': bytes.fromhex('F101'),
                'type': 'vertices'
            }
        }
        i = 0
        while i < len(self.dados) - 4:
            found_pattern = False
            for padrao_nome, padrao_info in padroes.items():
                prefixo = padrao_info['prefixo']
                if self.dados[i:i+2] == prefixo:
                    XX = self.dados[i+2]
                    if not (0x00 <= XX <= 0x2A):
                        continue  # XX fora do intervalo esperado
                    # Verifica o byte final do cabeçalho de vértices (deve ser 0x69)
                    if self.dados[i+3] != 0x69:
                        continue

                    found_pattern = True

                    grupo = {
                        'padrao': padrao_nome,
                        'inicio': i,
                        'XX': XX,
                        'vertices': [],
                        'faces_ativadas': [],
                        'faces_desativadas': [],
                        'uvs': [],  # Lista de pares (U,V)
                        'cor': self.gerar_cor_grupo()
                    }
                    # Lê os vértices
                    inicio_vertices = i + 4
                    quantidade_vertices_bytes = XX * 6  # Cada vértice ocupa 6 bytes
                    if inicio_vertices + quantidade_vertices_bytes > len(self.dados):
                        self.error_message = "Arquivo PCK corrompido ou incompleto (vértices)."
                        return False
                    for j in range(0, quantidade_vertices_bytes, 6):
                        try:
                            x, y, z = struct.unpack_from('<hhh', self.dados, inicio_vertices + j)
                            grupo['vertices'].append((x, y, z))
                        except struct.error as e:
                            self.error_message = f"Erro ao parsear vértices: {e}"
                            return False
                    grupo['hex_vertices'] = self.dados[inicio_vertices:inicio_vertices + quantidade_vertices_bytes].hex().upper()
                    grupo['vertices_offset'] = inicio_vertices
                    grupo['vertices_fim'] = inicio_vertices + quantidade_vertices_bytes

                    # Lê os UVs imediatamente após os vértices, se houver cabeçalho esperado
                    pos_after_vertices = grupo['vertices_fim']
                    if pos_after_vertices + 4 <= len(self.dados):
                        expected_uv_header = padrao_info['uv_prefixo'] + bytes([XX]) + bytes.fromhex('65')
                        if self.dados[pos_after_vertices: pos_after_vertices+4] == expected_uv_header:
                            inicio_uv = pos_after_vertices + 4
                            quantidade_uv_bytes = XX * 4  # Cada par UV ocupa 4 bytes
                            if inicio_uv + quantidade_uv_bytes > len(self.dados):
                                self.error_message = "Arquivo PCK corrompido ou incompleto (UVs)."
                                return False
                            for j in range(0, quantidade_uv_bytes, 4):
                                try:
                                    u, v = struct.unpack_from('<hh', self.dados, inicio_uv + j)
                                    grupo['uvs'].append((u, v))
                                except struct.error as e:
                                    self.error_message = f"Erro ao parsear UVs: {e}"
                                    return False
                            grupo['hex_uvs'] = self.dados[inicio_uv: inicio_uv + quantidade_uv_bytes].hex().upper()
                            grupo['uvs_offset'] = inicio_uv
                            grupo['uvs_fim'] = inicio_uv + quantidade_uv_bytes
                            pos_after_uv = grupo['uvs_fim']
                        else:
                            grupo['uvs'] = []
                            pos_after_uv = pos_after_vertices
                    else:
                        grupo['uvs'] = []
                        pos_after_uv = pos_after_vertices

                    # Procura o padrão de faces a partir de pos_after_uv
                    faces_prefixo = padrao_info['faces_prefixo'] + bytes([XX]) + bytes.fromhex('6A')
                    pos_faces = self.dados.find(faces_prefixo, pos_after_uv)
                    if pos_faces == -1:
                        self.error_message = "Padrão de faces não encontrado."
                        return False
                    inicio_faces = pos_faces + len(faces_prefixo) + 6  # Pula 6 bytes após o cabeçalho de faces
                    quantidade_faces_bytes = (XX - 2) * 3
                    if inicio_faces + quantidade_faces_bytes > len(self.dados):
                        self.error_message = "Arquivo PCK corrompido ou incompleto (faces)."
                        return False
                    for j in range(0, quantidade_faces_bytes, 3):
                        try:
                            face_bytes = struct.unpack_from('<BBB', self.dados, inicio_faces + j)
                            valor = face_bytes[0]
                            face_idx = j // 3 + 1
                            v1 = face_idx
                            v2 = face_idx + 1
                            v3 = face_idx + 2
                            if v3 > XX:
                                continue
                            face = {'indices': (v1, v2, v3), 'byte': valor}
                            if valor % 2 == 0:
                                grupo['faces_ativadas'].append(face)
                            else:
                                grupo['faces_desativadas'].append(face)
                        except struct.error as e:
                            self.error_message = f"Erro ao parsear faces: {e}"
                            return False
                    grupo['hex_faces'] = self.dados[pos_faces: pos_faces + len(faces_prefixo) + quantidade_faces_bytes].hex().upper()
                    grupo['faces_offset'] = pos_faces
                    grupo['faces_fim'] = pos_faces + len(faces_prefixo) + quantidade_faces_bytes

                    self.vertex_groups.append(grupo)
                    i = inicio_faces + quantidade_faces_bytes
                    break
            if not found_pattern:
                i += 1

        print(f"Total de Grupos de Vértices Encontrados: {len(self.vertex_groups)}")
        for idx, grupo in enumerate(self.vertex_groups):
            print(f"Grupo {idx + 1}:")
            print(f"  Padrão: {grupo['padrao']}")
            print(f"  Quantidade de Vértices: {grupo['XX']}")
            print(f"  Offset Vértices: {hex(grupo['vertices_offset'])} - {hex(grupo['vertices_fim'])}")
            print(f"  UVs Encontrados: {len(grupo['uvs'])}")
            print(f"  Faces Ativadas: {len(grupo['faces_ativadas'])}")
            print(f"  Faces Desativadas: {len(grupo['faces_desativadas'])}")

        total_vertices = sum(len(grupo['vertices']) for grupo in self.vertex_groups)
        total_faces_ativadas = sum(len(grupo['faces_ativadas']) for grupo in self.vertex_groups)
        total_faces_desativadas = sum(len(grupo['faces_desativadas']) for grupo in self.vertex_groups)
        print(f"Total de Vértices Globais: {total_vertices}")
        print(f"Total de Faces Ativadas Globais: {total_faces_ativadas}")
        print(f"Total de Faces Desativadas Globais: {total_faces_desativadas}")

        self.parse_success = True
        return True

    def gerar_cor_grupo(self):
        return (random.random(), random.random(), random.random(), 1)

    def parse(self):
        if not self.ler_arquivo_pck():
            return False
        sucesso = self.encontrar_padroes_vertices()
        if not sucesso:
            return False

        # Ajusta os índices globais de vértices para as faces
        self.faces_ativadas = []
        self.faces_desativadas = []
        self.vertices = []
        for grupo in self.vertex_groups:
            base_index = len(self.vertices)
            self.vertices.extend(grupo['vertices'])
            grupo['base_index'] = base_index
            grupo['faces_mesh_ativadas'] = []
            grupo['faces_mesh_desativadas'] = []
            for face in grupo['faces_ativadas']:
                if isinstance(face, dict) and 'indices' in face:
                    adjusted_indices = tuple(idx - 1 + base_index for idx in face['indices'])
                    if all(0 <= idx < len(self.vertices) for idx in adjusted_indices):
                        grupo['faces_mesh_ativadas'].append(adjusted_indices)
                else:
                    self.error_message = f"Formato de face inválido em grupo {grupo['padrao']}."
                    return False
            for face in grupo['faces_desativadas']:
                if isinstance(face, dict) and 'indices' in face:
                    adjusted_indices = tuple(idx - 1 + base_index for idx in face['indices'])
                    if all(0 <= idx < len(self.vertices) for idx in adjusted_indices):
                        grupo['faces_mesh_desativadas'].append(adjusted_indices)
                else:
                    self.error_message = f"Formato de face inválido em grupo {grupo['padrao']}."
                    return False

        print(f"Total de Vértices Globais: {len(self.vertices)}")
        total_faces_ativadas = sum(len(grupo['faces_mesh_ativadas']) for grupo in self.vertex_groups)
        total_faces_desativadas = sum(len(grupo['faces_mesh_desativadas']) for grupo in self.vertex_groups)
        print(f"Total de Faces Ativadas Globais: {total_faces_ativadas}")
        print(f"Total de Faces Desativadas Globais: {total_faces_desativadas}")

        self.parse_success = True
        return True


class InformationPanel(QDockWidget):
    def __init__(self, parent=None):
        super(InformationPanel, self).__init__("Informações", parent)
        self.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(['Propriedade', 'Valor'])
        self.setWidget(self.tree)

    def update_info(self, vertex_groups):
        self.tree.clear()
        for idx, grupo in enumerate(vertex_groups):
            grupo_item = QTreeWidgetItem(self.tree, [f'Grupo {idx + 1} - {grupo["padrao"]}', ''])
            cor = grupo['cor']
            qcolor = QColor(int(cor[0]*255), int(cor[1]*255), int(cor[2]*255), int(cor[3]*255))
            grupo_item.setBackground(0, qcolor)
            vertices_item = QTreeWidgetItem(grupo_item, ['Vértices', ''])
            QTreeWidgetItem(vertices_item, ['Quantidade', str(grupo['XX'])])
            QTreeWidgetItem(vertices_item, ['Tamanho (bytes)', str(grupo['vertices_fim'] - grupo['vertices_offset'])])
            QTreeWidgetItem(vertices_item, ['Offset Início', hex(grupo['vertices_offset'])])
            QTreeWidgetItem(vertices_item, ['Offset Fim', hex(grupo['vertices_fim'])])
            QTreeWidgetItem(vertices_item, ['Hex Vértices', grupo.get('hex_vertices', '')])
            uvs_item = QTreeWidgetItem(grupo_item, ['UVs', ''])
            QTreeWidgetItem(uvs_item, ['Quantidade', str(len(grupo['uvs']))])
            QTreeWidgetItem(uvs_item, ['Offset Início', hex(grupo.get('uvs_offset', 0))])
            QTreeWidgetItem(uvs_item, ['Offset Fim', hex(grupo.get('uvs_fim', 0))])
            QTreeWidgetItem(uvs_item, ['Hex UVs', grupo.get('hex_uvs', '')])
            faces_ativadas_item = QTreeWidgetItem(grupo_item, ['Faces Ativadas', ''])
            QTreeWidgetItem(faces_ativadas_item, ['Quantidade', str(len(grupo['faces_mesh_ativadas']))])
            for face_idx, face in enumerate(grupo['faces_mesh_ativadas']):
                face_info = f'Face {face_idx + 1}: {face}'
                QTreeWidgetItem(faces_ativadas_item, [face_info, 'Ativada'])
            faces_desativadas_item = QTreeWidgetItem(grupo_item, ['Faces Desativadas', ''])
            QTreeWidgetItem(faces_desativadas_item, ['Quantidade', str(len(grupo['faces_mesh_desativadas']))])
            for face_idx, face in enumerate(grupo['faces_mesh_desativadas']):
                face_info = f'Face {face_idx + 1}: {face}'
                QTreeWidgetItem(faces_desativadas_item, [face_info, 'Desativada'])
            total_faces = len(grupo['faces_mesh_ativadas']) + len(grupo['faces_mesh_desativadas'])
            QTreeWidgetItem(grupo_item, ['Faces Totais', str(total_faces)])
        self.tree.expandAll()


class TogglePanel(QDockWidget):
    def __init__(self, parent=None):
        super(TogglePanel, self).__init__("Controle de Faces Desativadas", parent)
        self.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
        widget = QWidget()
        layout = QVBoxLayout()
        self.checkbox = QCheckBox("Mostrar Faces Desativadas")
        self.checkbox.setChecked(False)
        self.checkbox.setToolTip("Marque para exibir as faces desativadas na visualização 3D.")
        layout.addWidget(self.checkbox)
        layout.addStretch()
        widget.setLayout(layout)
        self.setWidget(widget)


class SettingsPanel(QDockWidget):
    def __init__(self, parent=None):
        super(SettingsPanel, self).__init__("Configurações", parent)
        self.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
        widget = QWidget()
        layout = QVBoxLayout()
        bg_color_layout = QHBoxLayout()
        bg_color_label = QLabel("Cor de Fundo:")
        self.bg_color_button = QPushButton("Escolher Cor")
        self.bg_color_button.setToolTip("Clique para escolher a cor de fundo da visualização 3D.")
        bg_color_layout.addWidget(bg_color_label)
        bg_color_layout.addWidget(self.bg_color_button)
        layout.addLayout(bg_color_layout)
        self.wireframe_checkbox = QCheckBox("Mostrar Wireframe")
        self.wireframe_checkbox.setChecked(False)
        self.wireframe_checkbox.setToolTip("Ative para visualizar os meshes em wireframe.")
        layout.addWidget(self.wireframe_checkbox)
        self.grid_checkbox = QCheckBox("Mostrar Grade")
        self.grid_checkbox.setChecked(True)
        self.grid_checkbox.setToolTip("Ative para mostrar a grade de referência na cena 3D.")
        layout.addWidget(self.grid_checkbox)
        self.axes_checkbox = QCheckBox("Mostrar Eixos")
        self.axes_checkbox.setChecked(True)
        self.axes_checkbox.setToolTip("Ative para mostrar os eixos coordenados na cena 3D.")
        layout.addWidget(self.axes_checkbox)
        export_layout = QHBoxLayout()
        export_label = QLabel("Exportar Modelo:")
        self.export_button = QPushButton("Exportar para OBJ")
        self.export_button.setToolTip("Clique para exportar o modelo atual para um arquivo OBJ (incluindo UVs), omitindo faces desativadas.")
        export_layout.addWidget(export_label)
        export_layout.addWidget(self.export_button)
        layout.addLayout(export_layout)
        transparency_layout = QHBoxLayout()
        transparency_label = QLabel("Transparência do Mesh:")
        self.transparency_slider = QSlider(Qt.Horizontal)
        self.transparency_slider.setRange(0, 100)
        self.transparency_slider.setValue(100)
        self.transparency_slider.setToolTip("Ajuste a transparência dos meshes.")
        transparency_layout.addWidget(transparency_label)
        transparency_layout.addWidget(self.transparency_slider)
        layout.addLayout(transparency_layout)
        shading_layout = QHBoxLayout()
        shading_label = QLabel("Tipo de Shading:")
        self.shading_combo = QComboBox()
        self.shading_combo.addItems(["Flat", "Smooth"])
        self.shading_combo.setToolTip("Selecione o tipo de shading para os meshes.")
        shading_layout.addWidget(shading_label)
        shading_layout.addWidget(self.shading_combo)
        layout.addLayout(shading_layout)
        scale_layout = QHBoxLayout()
        scale_label = QLabel("Escalar Modelo:")
        self.scale_checkbox = QCheckBox("Dividir Vértices por 256")
        self.scale_checkbox.setToolTip("Ative para escalar o modelo dividindo todas as coordenadas dos vértices por 256.")
        scale_layout.addWidget(scale_label)
        scale_layout.addWidget(self.scale_checkbox)
        layout.addLayout(scale_layout)
        layout.addStretch()
        widget.setLayout(layout)
        self.setWidget(widget)


class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.setWindowTitle("Visualizador 3D de Arquivos PCK")
        self.setGeometry(100, 100, 1200, 800)
        self.view_widget = QtInteractor(self)
        self.setCentralWidget(self.view_widget)

        # Melhora a visualização 3D: ativa anti-aliasing e configura o render window
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
        self.toggle_panel.checkbox.stateChanged.connect(self.toggle_faces_desativadas)
        self.settings_panel.bg_color_button.clicked.connect(self.choose_background_color)
        self.settings_panel.wireframe_checkbox.stateChanged.connect(self.toggle_wireframe)
        self.settings_panel.grid_checkbox.stateChanged.connect(self.toggle_grid)
        self.settings_panel.axes_checkbox.stateChanged.connect(self.toggle_axes)
        self.settings_panel.export_button.clicked.connect(self.export_to_obj)
        self.settings_panel.transparency_slider.valueChanged.connect(self.adjust_transparency)
        self.settings_panel.shading_combo.currentIndexChanged.connect(self.change_shading)
        self.settings_panel.scale_checkbox.stateChanged.connect(self.toggle_scale)
        self.mesh_grupos = []
        self.mesh_desativados = None
        self.current_shading = "Smooth"
        self.scale_enabled = False
        self.original_vertices = []
        self.init_menu()
        self.parser = None

    def init_menu(self):
        menubar = self.menuBar()
        file_menu = menubar.addMenu('Arquivo')
        abrir_action = QAction('Abrir PCK', self)
        abrir_action.setShortcut('Ctrl+O')
        abrir_action.setStatusTip('Abrir arquivo PCK')
        abrir_action.triggered.connect(self.abrir_arquivo)
        file_menu.addAction(abrir_action)
        exportar_action = QAction('Exportar para OBJ', self)
        exportar_action.setShortcut('Ctrl+E')
        exportar_action.setStatusTip('Exportar modelo para OBJ')
        exportar_action.triggered.connect(self.export_to_obj)
        file_menu.addAction(exportar_action)
        sair_action = QAction('Sair', self)
        sair_action.setShortcut('Ctrl+Q')
        sair_action.setStatusTip('Sair da aplicação')
        sair_action.triggered.connect(self.close)
        file_menu.addAction(sair_action)

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
        grid_actor = self.view_widget.add_mesh(grid, color='gray', line_width=1, opacity=0.5, name="Grade")
        return grid_actor

    def add_axes(self):
        axes_actor = self.view_widget.add_axes()
        return axes_actor

    def abrir_arquivo(self):
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly
        filepath, _ = QFileDialog.getOpenFileName(
            self,
            "Abrir Arquivo PCK",
            "",
            "PCK Files (*.pck);;All Files (*)",
            options=options
        )
        if filepath:
            parser = PCKParser(filepath)
            sucesso = parser.parse()
            if sucesso:
                self.original_vertices = parser.vertices.copy()
                self.render_model(parser.vertex_groups, parser.vertices)
                self.info_panel.update_info(parser.vertex_groups)
                self.setWindowTitle(f"Visualizador 3D de Arquivos PCK - {filepath}")
                self.parser = parser
            else:
                QMessageBox.critical(self, "Erro", parser.error_message)

    def render_model(self, vertex_groups, vertices):
        self.view_widget.clear()
        self.grid_actor = self.add_grid()
        self.axes_actor = self.add_axes()
        if self.scale_enabled:
            vertices_np = np.array(vertices, dtype=np.float32) / 256.0
        else:
            vertices_np = np.array(vertices, dtype=np.float32)
        self.mesh_grupos = []
        for idx, grupo in enumerate(vertex_groups):
            if grupo['faces_mesh_ativadas']:
                faces_np = np.array(grupo['faces_mesh_ativadas'], dtype=np.int32)
                faces_formatted = np.hstack([np.full((faces_np.shape[0], 1), 3), faces_np]).flatten()
                mesh = pv.PolyData(vertices_np, faces_formatted)
                mesh_plot = self.view_widget.add_mesh(
                    mesh,
                    color=grupo['cor'][:3],
                    show_edges=self.settings_panel.wireframe_checkbox.isChecked(),
                    edge_color='black',
                    name=f"Grupo {idx + 1}",
                    smooth_shading=(self.current_shading == "Smooth")
                )
                self.mesh_grupos.append(mesh_plot)
        self.mesh_desativados = None
        faces_desativadas_all = []
        for grupo in vertex_groups:
            if grupo['faces_mesh_desativadas']:
                faces_desativadas_all.extend(grupo['faces_mesh_desativadas'])
        if faces_desativadas_all:
            faces_desativadas_np = np.array(faces_desativadas_all, dtype=np.int32)
            faces_formatted = np.hstack([np.full((faces_desativadas_np.shape[0], 1), 3), faces_desativadas_np]).flatten()
            mesh_desativado = pv.PolyData(vertices_np, faces_formatted)
            self.mesh_desativados = self.view_widget.add_mesh(
                mesh_desativado,
                color=(1, 0, 0),
                show_edges=self.settings_panel.wireframe_checkbox.isChecked(),
                edge_color='black',
                name="Faces Desativadas",
                smooth_shading=(self.current_shading == "Smooth")
            )
            self.mesh_desativados.SetVisibility(False)
        self.view_widget.reset_camera()

    def toggle_faces_desativadas(self, state):
        if self.mesh_desativados:
            self.mesh_desativados.SetVisibility(state == Qt.Checked)
            self.view_widget.render()

    def choose_background_color(self):
        color = QColorDialog.getColor()
        if color.isValid():
            rgb_color = (color.red() / 255, color.green() / 255, color.blue() / 255)
            self.view_widget.set_background(rgb_color)

    def toggle_wireframe(self, state):
        show_wireframe = state == Qt.Checked
        for mesh in self.mesh_grupos:
            mesh.GetProperty().SetEdgeVisibility(show_wireframe)
        if self.mesh_desativados:
            self.mesh_desativados.GetProperty().SetEdgeVisibility(show_wireframe)
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
        for mesh in self.mesh_grupos:
            # Ajusta a opacidade usando a propriedade do ator VTK
            mesh.GetProperty().SetOpacity(transparency)
        if self.mesh_desativados:
            self.mesh_desativados.GetProperty().SetOpacity(transparency)
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
            QMessageBox.warning(self, "Aviso", "Nenhum modelo carregado para exportação.")
            return

        filepath, _ = QFileDialog.getSaveFileName(
            self,
            "Exportar para OBJ",
            "",
            "OBJ Files (*.obj);;All Files (*)"
        )
        if not filepath:
            return

        if self.scale_enabled:
            vertices = np.array(self.original_vertices, dtype=np.float32) / 256.0
        else:
            vertices = np.array(self.original_vertices, dtype=np.float32)

        # Constrói os UVs agrupados de cada grupo, na ordem dos vértices
        global_uvs = []
        for grupo in self.parser.vertex_groups:
            if len(grupo['uvs']) == len(grupo['vertices']):
                global_uvs.extend(grupo['uvs'])
            else:
                global_uvs.extend([(0, 0)] * len(grupo['vertices']))
        if len(global_uvs) < len(vertices):
            global_uvs.extend([(0, 0)] * (len(vertices) - len(global_uvs)))
        elif len(global_uvs) > len(vertices):
            global_uvs = global_uvs[:len(vertices)]

        obj_lines = []
        for v in vertices:
            obj_lines.append(f"v {v[0]} {v[1]} {v[2]}")
        for uv in global_uvs:
            obj_lines.append(f"vt {uv[0]} {uv[1]}")
        for grupo in self.parser.vertex_groups:
            for face in grupo['faces_mesh_ativadas']:
                v1, v2, v3 = face
                obj_lines.append(f"f {v1+1}/{v1+1} {v2+1}/{v2+1} {v3+1}/{v3+1}")
        obj_content = "\n".join(obj_lines)
        print("Conteúdo do arquivo OBJ:")
        print(obj_content)
        try:
            with open(filepath, "w") as f:
                f.write(obj_content)
            QMessageBox.information(self, "Sucesso", f"Modelo exportado com sucesso para {filepath}.")
        except Exception as e:
            QMessageBox.critical(self, "Erro", f"Falha ao exportar o modelo: {e}")

    def change_lighting(self, index):
        pass

def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
