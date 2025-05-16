import sys
import os
import pandas as pd
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QTabWidget, QLabel, QLineEdit, 
                             QPushButton, QTextEdit, QGroupBox, QGridLayout, 
                             QComboBox, QFileDialog, QCheckBox, QProgressBar,
                             QListWidget, QMessageBox, QAction, QToolBar)
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtGui import QIcon
from datetime import datetime, timedelta
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.io.sac import SACTrace

# 导入项目中的模块
sys.path.append(os.path.join(os.path.dirname(__file__), 'Python'))
import a1_get_metadata
import a2_download_daydata
import a3_get_catalog
import a4_download_event
import a4_cut_event

# 导入进度条适配器
from progress_adapter import TqdmProgressAdapter

class WorkerThread(QThread):
    """工作线程，用于执行耗时的下载任务"""
    progressUpdate = pyqtSignal(int)
    statusUpdate = pyqtSignal(str)
    finished = pyqtSignal()
    
    def __init__(self, task_type, params):
        super().__init__()
        self.task_type = task_type
        self.params = params
        
    def run(self):
        try:
            if self.task_type == "metadata":
                self.statusUpdate.emit("正在获取台站元数据...")
                # 设置全局变量
                a1_get_metadata.network_code = self.params["network_code"]
                a1_get_metadata.network_starttime = self.params["network_starttime"]
                a1_get_metadata.network_endtime = self.params["network_endtime"]
                a1_get_metadata.channel_query = self.params["channel_query"]
                
                # 重新获取inventory
                self.statusUpdate.emit(f"连接IRIS服务器...")
                client = a1_get_metadata.client
                
                self.statusUpdate.emit(f"获取网络 {a1_get_metadata.network_code} 的台站信息...")
                inventory = client.get_stations(network=a1_get_metadata.network_code,
                                                channel=a1_get_metadata.channel_query,
                                                starttime=a1_get_metadata.network_starttime,
                                                endtime=a1_get_metadata.network_endtime,
                                                level='channel')
                
                # 处理inventory
                for network in inventory:
                    # 重置变量
                    a1_get_metadata.stations_code = []
                    a1_get_metadata.stations_starttime = []
                    a1_get_metadata.stations_endtime = []
                    a1_get_metadata.stations_channels = []
                    a1_get_metadata.stations_total_channels = []
                    a1_get_metadata.stations_site = []
                    a1_get_metadata.stations_latitude = []
                    a1_get_metadata.stations_longitude = []
                    a1_get_metadata.stations_elevation = []
                    
                    # 更新进度条最大值
                    total_stations = len(network.stations)
                    self.statusUpdate.emit(f"找到 {total_stations} 个台站")
                    
                    # 处理每个台站
                    for i, station in enumerate(network.stations):
                        self.statusUpdate.emit(f"处理台站: {station.code} ({i+1}/{total_stations})")
                        a1_get_metadata.process_station(station)
                        # 更新进度
                        progress = int(100 * (i + 1) / total_stations)
                        self.progressUpdate.emit(progress)
                    
                    # 创建DataFrame并保存
                    df = a1_get_metadata.create_dataframe()
                    filename = f'{a1_get_metadata.network_code}_fdsn_metadata.txt'
                    df.to_csv(filename, sep='\t', float_format='%.6f')
                    self.statusUpdate.emit(f"元数据已保存到 {filename}")
                
            elif self.task_type == "daydata":
                self.statusUpdate.emit("正在下载日数据...")
                # 设置参数
                output_dir = self.params["output_dir"]
                network = self.params["network"]
                location = self.params["location"]
                stations_download = self.params["stations_download"]
                output_units_seis = self.params["output_units_seis"]
                output_units_pres = self.params["output_units_pres"]
                metadatafile = self.params["metadatafile"]
                
                # 创建目录
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                
                if not os.path.exists(f'{output_dir}/{network}'):
                    os.makedirs(f'{output_dir}/{network}')

                # 读取元数据
                metadata = pd.read_csv(metadatafile, sep='\t', index_col=0)
                
                stations = metadata['station']
                stations_starttime = metadata['starttime']
                stations_endtime = metadata['endtime']
                stations_latitude = metadata['latitude']
                stations_longitude = metadata['longitude']
                stations_elevation = metadata['elevation']
                stations_channels = metadata['channel']
                
                # 连接IRIS服务器
                self.statusUpdate.emit("连接IRIS服务器...")
                client = Client('IRIS')
                
                # 计算总数据量
                data_numbers = 0
                for ista, station in enumerate(stations):
                    if station not in stations_download:
                        continue
                    
                    channels_str = stations_channels[ista]
                    channels = channels_str.split(',')
                    
                    start_date = datetime.strptime(stations_starttime[ista][0:10],'%Y-%m-%d')
                    end_date = datetime.strptime(stations_endtime[ista][0:10],'%Y-%m-%d')
                    
                    delta = abs(end_date - start_date)
                    data_number = len(channels) * delta.days
                    data_numbers += data_number
                
                self.statusUpdate.emit(f"总共需要下载 {data_numbers} 个数据文件")
                
                # 创建进度条适配器
                progress_adapter = TqdmProgressAdapter(
                    signal_progress=self.progressUpdate,
                    signal_status=self.statusUpdate,
                    total=data_numbers
                )
                
                # 开始下载
                current_progress = 0
                for ista, station in enumerate(stations):
                    if station not in stations_download:
                        continue
                    
                    if not os.path.exists(f'{output_dir}/{network}/{station}'):
                        os.makedirs(f'{output_dir}/{network}/{station}')
                    
                    sta_lat = stations_latitude[ista]
                    sta_lon = stations_longitude[ista]
                    sta_ele = stations_elevation[ista]
                    channels_str = stations_channels[ista]
                    
                    channels = channels_str.split(',')
                    
                    start_date = datetime.strptime(stations_starttime[ista][0:10],'%Y-%m-%d')
                    end_date = datetime.strptime(stations_endtime[ista][0:10],'%Y-%m-%d')
                    idate = start_date
                    
                    while idate <= end_date:
                        dayid = datetime.strftime(idate,'%Y%m%d_%H%M%S')
                        
                        starttime_str = datetime.strftime(idate,'%Y-%m-%dT%H:%M:%S')
                        endtime_str = datetime.strftime(idate+timedelta(days=1),'%Y-%m-%dT%H:%M:%S')
                        starttime = UTCDateTime(starttime_str)
                        endtime = UTCDateTime(endtime_str)
                        
                        idate = idate + timedelta(days=1)
                        
                        for channel in channels:
                            # 更新进度
                            progress_adapter.update(1)
                            current_progress += 1
                            
                            filename = f'{output_dir}/{network}/{station}/{dayid}_{network}_{station}_{channel}.SAC'
                            
                            if os.path.isfile(filename):
                                self.statusUpdate.emit(f'{filename} 已存在. 跳过!')
                                continue
                            
                            # 下载波形数据
                            self.statusUpdate.emit(f'下载台站: {station} 通道: {channel} 从: {starttime_str} 到: {endtime_str}')
                            try:
                                st = client.get_waveforms(network=network, station=station, location=location, channel=channel,
                                                          starttime=starttime, endtime=endtime, attach_response=True)
                                
                                # 预处理
                                if channel[1] == 'H':  # 地震计
                                    st.remove_response(output=output_units_seis)
                                elif channel[-1] == 'H':  # DPG或水听器
                                    st.remove_response(output=output_units_pres)
                                
                                # 去均值
                                st.detrend('demean')
                                
                                # 写入SAC文件
                                st.write(filename, format='SAC')
                                
                                # 读取仅头部
                                sac = SACTrace.read(filename, headonly=True)
                                
                                # 写入SAC头部信息
                                sac.knetwk = network
                                sac.kstnm = station
                                sac.khole = location
                                sac.kcmpnm = channel
                                
                                sac.stla = sta_lat
                                sac.stlo = sta_lon
                                sac.stel = sta_ele
                                
                                # 只写入头部，文件必须存在
                                sac.write(filename, headonly=True)
                                
                                self.statusUpdate.emit(f'文件 {filename} 下载并保存成功')
                                
                            except Exception as e:
                                self.statusUpdate.emit(f'错误: {str(e)}')
                                self.statusUpdate.emit(f'无法下载 {filename}. 跳过!')
                
                # 关闭进度条
                progress_adapter.close()
                self.statusUpdate.emit("日数据下载完成！")
                
            elif self.task_type == "catalog":
                self.statusUpdate.emit("正在获取地震目录...")
                # 设置参数
                output_dir = self.params["output_dir"]
                network = self.params["network"]
                minmag_local = self.params["minmag_local"]
                maxmag_local = self.params["maxmag_local"]
                maxradius_local = self.params["maxradius_local"]
                minmag_global = self.params["minmag_global"]
                metadatafile = self.params["metadatafile"]
                
                # 创建目录
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                
                if not os.path.exists(f'{output_dir}/{network}'):
                    os.makedirs(f'{output_dir}/{network}')
                
                # 读取元数据
                metadata = pd.read_csv(metadatafile, sep='\t', index_col=0)
                
                stations = metadata['station']
                stations_starttime = metadata['starttime']
                stations_endtime = metadata['endtime']
                stations_latitude = metadata['latitude']
                stations_longitude = metadata['longitude']
                
                # 连接IRIS服务器
                self.statusUpdate.emit("连接IRIS服务器...")
                client = Client('IRIS')
                
                # 创建进度条
                total_stations = len(stations)
                self.statusUpdate.emit(f"总共需要获取 {total_stations} 个台站的地震目录")
                
                # 开始获取地震目录
                for ista in range(len(stations)):
                    station = stations[ista]
                    self.statusUpdate.emit(f"正在获取台站 {station} 的地震目录 ({ista+1}/{total_stations})")
                    
                    sta_lat = stations_latitude[ista]
                    sta_lon = stations_longitude[ista]
                    
                    starttime = UTCDateTime(stations_starttime[ista])
                    endtime = UTCDateTime(stations_endtime[ista])
                    
                    try:
                        # 获取本地地震目录
                        self.statusUpdate.emit(f"获取台站 {station} 周围的本地地震目录 (震级 {minmag_local}-{maxmag_local}，半径 {maxradius_local}°)")
                        cata_local = client.get_events(
                            starttime=starttime, 
                            endtime=endtime, 
                            minmagnitude=minmag_local, 
                            maxmagnitude=maxmag_local,
                            latitude=sta_lat, 
                            longitude=sta_lon, 
                            maxradius=maxradius_local
                        )
                        
                        # 获取全球地震目录
                        self.statusUpdate.emit(f"获取全球地震目录 (震级 >= {minmag_global})")
                        cata_global = client.get_events(
                            starttime=starttime, 
                            endtime=endtime, 
                            minmag=minmag_global
                        )
                        
                        # 合并地震目录
                        catalog = cata_local + cata_global
                        self.statusUpdate.emit(f"找到 {len(catalog)} 个地震事件")
                        
                        eventids = []
                        otimes_str = []
                        magnitude_types = []
                        descriptions = []
                        events_latitude = []
                        events_longitude = []
                        events_depth = []
                        magnitudes = []
                        
                        # 处理每个地震
                        for event in catalog:
                            eventid = event.resource_id.id.split('eventid=')[-1]
                            origin = event.origins[0]
                            otime_str = origin.time.__unicode__()
                            event_latitude = origin.latitude
                            event_longitude = origin.longitude
                            event_depth = origin.depth * 0.001  # 转换为公里
                            magnitude_type = event.magnitudes[0].magnitude_type
                            magnitude = event.magnitudes[0].mag
                            description = event.event_descriptions[0].text
                            
                            eventids.append(eventid)
                            otimes_str.append(otime_str)
                            magnitude_types.append(magnitude_type)
                            descriptions.append(description)
                            events_latitude.append(event_latitude)
                            events_longitude.append(event_longitude)
                            events_depth.append(event_depth)
                            magnitudes.append(magnitude)
                        
                        # 创建DataFrame并保存
                        df = pd.DataFrame({
                            'event id': eventids, 
                            'origin time': otimes_str, 
                            'latitude': events_latitude, 
                            'longitude': events_longitude, 
                            'depth': events_depth, 
                            'magnitude type': magnitude_types, 
                            'magnitude': magnitudes, 
                            'description': descriptions
                        })
                        
                        filename = f'{output_dir}/{network}/{network}_{station}_catalog.txt'
                        df.to_csv(filename, sep='\t', float_format='%.6f')
                        self.statusUpdate.emit(f"台站 {station} 的地震目录已保存到 {filename}")
                        
                    except Exception as e:
                        self.statusUpdate.emit(f"错误: {str(e)}")
                        self.statusUpdate.emit(f"无法获取台站 {station} 的地震目录. 跳过!")
                    
                    # 更新进度
                    progress = int(100 * (ista + 1) / total_stations)
                    self.progressUpdate.emit(progress)
                
                self.statusUpdate.emit("地震目录获取完成！")
                
            elif self.task_type == "eventdata":
                self.statusUpdate.emit("正在下载事件数据...")
                # 设置参数
                output_dir = self.params["output_dir"]
                network = self.params["network"]
                location = self.params["location"]
                stations_download = self.params["stations_download"]
                output_units_seis = self.params["output_units_seis"]
                output_units_pres = self.params["output_units_pres"]
                event_length = self.params["event_length"]
                btime = self.params["btime"]
                Rayleigh_velocity = self.params["Rayleigh_velocity"]
                metadatafile = self.params["metadatafile"]
                
                # 创建目录
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                
                # 读取元数据
                metadata = pd.read_csv(metadatafile, sep='\t', index_col=0)
                
                stations = metadata['station']
                stations_starttime = metadata['starttime']
                stations_endtime = metadata['endtime']
                stations_latitude = metadata['latitude']
                stations_longitude = metadata['longitude']
                stations_elevation = metadata['elevation']
                stations_channels = metadata['channel']
                
                # 初始化TauPy模型和IRIS客户端
                from obspy.taup import TauPyModel
                from obspy.geodetics import gps2dist_azimuth, kilometer2degrees
                
                self.statusUpdate.emit("初始化TauPy模型...")
                model = TauPyModel(model='iasp91')
                
                self.statusUpdate.emit("连接IRIS服务器...")
                client = Client('IRIS')
                
                # 计算需要下载的数据总量
                data_numbers = 0
                for ista, station in enumerate(stations):
                    if station not in stations_download:
                        continue
                    
                    channels_str = stations_channels[ista]        
                    channels = channels_str.split(',')
                    
                    # 读取台站的地震目录
                    catalogfile = f'catalog/{network}/{network}_{station}_catalog.txt'
                    if not os.path.exists(catalogfile):
                        self.statusUpdate.emit(f"警告: 台站 {station} 的地震目录文件不存在: {catalogfile}")
                        continue
                        
                    catalog = pd.read_csv(catalogfile, sep='\t', index_col=0)
                    
                    data_number = len(channels) * len(catalog)
                    data_numbers += data_number
                
                self.statusUpdate.emit(f"总共需要下载 {data_numbers} 个事件数据文件")
                
                # 创建进度条适配器
                progress_adapter = TqdmProgressAdapter(
                    signal_progress=self.progressUpdate,
                    signal_status=self.statusUpdate,
                    total=data_numbers
                )
                
                # 开始下载事件数据
                for ista, station in enumerate(stations):
                    if station not in stations_download:
                        continue
                    
                    sta_lat = stations_latitude[ista]
                    sta_lon = stations_longitude[ista]
                    sta_ele = stations_elevation[ista]
                    channels_str = stations_channels[ista]
                    
                    channels = channels_str.split(',')
                    
                    # 读取台站的地震目录
                    catalogfile = f'catalog/{network}/{network}_{station}_catalog.txt'
                    if not os.path.exists(catalogfile):
                        continue
                        
                    catalog = pd.read_csv(catalogfile, sep='\t', index_col=0)
                    
                    otimes_str = catalog['origin time']
                    events_latitude = catalog['latitude']
                    events_longitude = catalog['longitude']
                    events_depth = catalog['depth']
                    magnitude_types = catalog['magnitude type']
                    magnitudes = catalog['magnitude']
                    
                    for ievt, otime_str in enumerate(otimes_str):
                        otime = UTCDateTime(otime_str)
                        evt_lat = events_latitude[ievt]
                        evt_lon = events_longitude[ievt]
                        evt_dep = events_depth[ievt]
                        magnitude = magnitudes[ievt]
                        magnitude_type = magnitude_types[ievt]
                        
                        eventid = otime.datetime.strftime('%Y%m%d_%H%M%S')
                        
                        if not os.path.exists(f'{output_dir}/{eventid}'):
                            os.makedirs(f'{output_dir}/{eventid}')
                            
                        for channel in channels:
                            # 更新进度
                            progress_adapter.update(1)
                            
                            filename = f'{output_dir}/{eventid}/{eventid}_{network}_{station}_{channel}.SAC'
                            
                            if os.path.isfile(filename):
                                self.statusUpdate.emit(f'{filename} 已存在. 跳过!')
                                continue
                                
                            starttime = otime + btime
                            endtime = otime + event_length + btime
                            
                            starttime_str = starttime.datetime.strftime('%Y-%m-%dT%H:%M:%S')
                            endtime_str = endtime.datetime.strftime('%Y-%m-%dT%H:%M:%S')
                            
                            # 下载波形数据
                            self.statusUpdate.emit(f'下载台站: {station} 通道: {channel} 从: {starttime_str} 到: {endtime_str}')
                            try:
                                st = client.get_waveforms(network=network, station=station, location=location, channel=channel,
                                                          starttime=starttime, endtime=endtime, attach_response=True)
                                
                                # 预处理
                                if channel[1] == 'H':  # 地震计
                                    st.remove_response(output=output_units_seis)
                                elif channel[-1] == 'H':  # DPG或水听器
                                    st.remove_response(output=output_units_pres)
                                
                                # 去均值和锥形窗处理
                                st.detrend('demean')
                                st.taper(max_percentage=0.1)
                                
                                # 写入SAC文件
                                st.write(filename, format='SAC')
                                
                                # 计算距离和方位角
                                distance_in_m, baz, az = gps2dist_azimuth(sta_lat, sta_lon, evt_lat, evt_lon)
                                distance_in_km = distance_in_m * 0.001
                                distance_in_degree = kilometer2degrees(distance_in_km)
                                
                                # 获取P波和S波到时
                                P_arrivals = model.get_travel_times(source_depth_in_km=evt_dep, 
                                                                    distance_in_degree=distance_in_degree, 
                                                                    phase_list='P')
                                S_arrivals = model.get_travel_times(source_depth_in_km=evt_dep, 
                                                                    distance_in_degree=distance_in_degree, 
                                                                    phase_list='S')
                                
                                # 瑞利波窗口
                                Rayleigh_begin = distance_in_km / Rayleigh_velocity[1]
                                Rayleigh_end = distance_in_km / Rayleigh_velocity[0]
                                
                                # 读取仅头部
                                sac = SACTrace.read(filename, headonly=True)
                                
                                # 写入SAC头部信息
                                sac.knetwk = network
                                sac.kstnm = station
                                sac.khole = location
                                sac.kcmpnm = channel
                                
                                sac.stla = sta_lat
                                sac.stlo = sta_lon
                                sac.stel = sta_ele
                                
                                sac.evla = evt_lat
                                sac.evlo = evt_lon
                                sac.evdp = evt_dep
                                sac.mag = magnitude
                                
                                sac.gcarc = distance_in_degree
                                sac.dist = distance_in_km
                                sac.az = az
                                sac.baz = baz
                                
                                sac.o = otime
                                sac.iztype = 'io'
                                sac.ko = 'O'
                                
                                # P波和S波到时
                                if P_arrivals:
                                    sac.a = otime + P_arrivals[0].time
                                    sac.t0 = otime + S_arrivals[0].time
                                    sac.ka = 'P'
                                    sac.kt0 = 'S'
                                
                                # 瑞利波窗口
                                sac.t1 = otime + Rayleigh_begin
                                sac.t2 = otime + Rayleigh_end
                                sac.kt1 = 'RayStart'
                                sac.kt2 = 'RayEnd'
                                
                                # 写入发震时刻到SAC头部
                                sac.kevnm = otime_str[0:4] + otime_str[5:7] + otime_str[8:10] + \
                                            otime_str[10:13] + otime_str[14:16] + otime_str[17:19]
                                # 发震时刻的小数部分
                                sac.kuser0 = otime_str[19:27]
                                # 震级类型
                                sac.kuser1 = magnitude_type
                                
                                # 只写入头部，文件必须存在
                                sac.write(filename, headonly=True)
                                
                                self.statusUpdate.emit(f'文件 {filename} 下载并保存成功')
                                
                            except Exception as e:
                                self.statusUpdate.emit(f'错误: {str(e)}')
                                self.statusUpdate.emit(f'无法下载 {filename}. 跳过!')
                
                # 关闭进度条
                progress_adapter.close()
                self.statusUpdate.emit("事件数据下载完成！")
                
            elif self.task_type == "cutevents":
                self.statusUpdate.emit("正在从日数据中截取事件数据...")
                # 设置参数
                input_dir = self.params["input_dir"]
                output_dir = self.params["output_dir"]
                catlog_dir = self.params["catlog_dir"]
                network = self.params["network"]
                location = self.params["location"]
                stations_download = self.params["stations_download"]
                event_length = self.params["event_length"]
                btime = self.params["btime"]
                Rayleigh_velocity = self.params["Rayleigh_velocity"]
                metadatafile = self.params["metadatafile"]
                
                # 创建目录
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                
                # 读取元数据
                metadata = pd.read_csv(metadatafile, sep='\t', index_col=0)
                
                stations = metadata['station']
                stations_starttime = metadata['starttime']
                stations_endtime = metadata['endtime']
                stations_latitude = metadata['latitude']
                stations_longitude = metadata['longitude']
                stations_elevation = metadata['elevation']
                stations_channels = metadata['channel']
                
                # 初始化TauPy模型
                from obspy.taup import TauPyModel
                from obspy.geodetics import gps2dist_azimuth, kilometer2degrees
                import glob
                from obspy import read
                
                self.statusUpdate.emit("初始化TauPy模型...")
                model = TauPyModel(model='iasp91')
                
                # 计算需要截取的数据总量
                data_numbers = 0
                for ista, station in enumerate(stations):
                    if station not in stations_download:
                        continue
                    
                    channels_str = stations_channels[ista]        
                    channels = channels_str.split(',')
                    
                    # 读取台站的地震目录
                    catalogfile = f'{catlog_dir}/{network}/{network}_{station}_catalog.txt'
                    if not os.path.exists(catalogfile):
                        self.statusUpdate.emit(f"警告: 台站 {station} 的地震目录文件不存在: {catalogfile}")
                        continue
                        
                    catalog = pd.read_csv(catalogfile, sep='\t', index_col=0)
                    
                    data_number = len(channels) * len(catalog)
                    data_numbers += data_number
                
                self.statusUpdate.emit(f"总共需要截取 {data_numbers} 个事件数据文件")
                
                # 创建进度条适配器
                progress_adapter = TqdmProgressAdapter(
                    signal_progress=self.progressUpdate,
                    signal_status=self.statusUpdate,
                    total=data_numbers
                )
                
                # 开始截取事件数据
                for ista, station in enumerate(stations):
                    if station not in stations_download:
                        continue
                    
                    sta_lat = stations_latitude[ista]
                    sta_lon = stations_longitude[ista]
                    sta_ele = stations_elevation[ista]
                    channels_str = stations_channels[ista]
                    
                    channels = channels_str.split(',')
                    
                    # 读取台站的地震目录
                    catalogfile = f'{catlog_dir}/{network}/{network}_{station}_catalog.txt'
                    if not os.path.exists(catalogfile):
                        continue
                        
                    catalog = pd.read_csv(catalogfile, sep='\t', index_col=0)
                    
                    otimes_str = catalog['origin time']
                    events_latitude = catalog['latitude']
                    events_longitude = catalog['longitude']
                    events_depth = catalog['depth']
                    magnitude_types = catalog['magnitude type']
                    magnitudes = catalog['magnitude']
                    
                    for ievt, otime_str in enumerate(otimes_str):
                        otime = UTCDateTime(otime_str)
                        evt_lat = events_latitude[ievt]
                        evt_lon = events_longitude[ievt]
                        evt_dep = events_depth[ievt]
                        magnitude = magnitudes[ievt]
                        magnitude_type = magnitude_types[ievt]
                        
                        eventid = otime.datetime.strftime('%Y%m%d_%H%M%S')
                        
                        if not os.path.exists(f'{output_dir}/{eventid}'):
                            os.makedirs(f'{output_dir}/{eventid}')
                            
                        starttime = otime + btime
                        endtime = otime + event_length + btime
                        
                        start_date = starttime.datetime.strftime('%Y%m%d')
                        end_date = endtime.datetime.strftime('%Y%m%d')
                        
                        starttime_str = starttime.datetime.strftime('%Y-%m-%dT%H:%M:%S')
                        endtime_str = endtime.datetime.strftime('%Y-%m-%dT%H:%M:%S')
                            
                        for channel in channels:
                            # 更新进度
                            progress_adapter.update(1)
                            
                            self.statusUpdate.emit(f'截取台站: {station} 通道: {channel} 从: {starttime_str} 到: {endtime_str}')
                            
                            filename = f'{output_dir}/{eventid}/{eventid}_{network}_{station}_{channel}.SAC'
                            
                            if os.path.isfile(filename):
                                self.statusUpdate.emit(f'{filename} 已存在. 跳过!')
                                continue
                            
                            # 检查是否有对应的日数据
                            if start_date == end_date:
                                day_filename = glob.glob(f'{input_dir}/{network}/{station}/*{start_date}*{channel}.SAC')
                                
                                if not day_filename:
                                    self.statusUpdate.emit('日数据不存在. 跳过!')
                                    continue
                                    
                                st = read(day_filename[0])
                            else:
                                start_filename = glob.glob(f'{input_dir}/{network}/{station}/*{start_date}*{channel}.SAC')
                                end_filename = glob.glob(f'{input_dir}/{network}/{station}/*{end_date}*{channel}.SAC')
                                
                                if not start_filename or not end_filename:
                                    self.statusUpdate.emit('日数据不存在. 跳过!')
                                    continue
                                    
                                st = read(start_filename[0])
                                st += read(end_filename[0])
                                st.merge()
                            
                            # 确保数据包含了整个事件
                            if st[0].stats.starttime > starttime or st[0].stats.endtime < endtime:
                                self.statusUpdate.emit('日数据不包含完整事件数据. 跳过!')
                                continue
                            
                            # 处理可能的掩码数据
                            if hasattr(st[0].data, 'filled') and callable(getattr(st[0].data, 'filled')):
                                st[0].data = st[0].data.filled()
                            
                            tr = st[0]
                            dt = tr.stats.delta
                            st_event = st.trim(starttime=starttime, endtime=endtime-dt)
                            
                            # 预处理
                            st_event.detrend('demean')
                            st_event.taper(max_percentage=0.1)
                            
                            # 写入SAC文件
                            st_event.write(filename, format='SAC')
                            
                            # 计算距离和方位角
                            distance_in_m, baz, az = gps2dist_azimuth(sta_lat, sta_lon, evt_lat, evt_lon)
                            distance_in_km = distance_in_m * 0.001
                            distance_in_degree = kilometer2degrees(distance_in_km)
                            
                            # 获取P波和S波到时
                            P_arrivals = model.get_travel_times(source_depth_in_km=evt_dep, 
                                                                distance_in_degree=distance_in_degree, 
                                                                phase_list='P')
                            S_arrivals = model.get_travel_times(source_depth_in_km=evt_dep, 
                                                                distance_in_degree=distance_in_degree, 
                                                                phase_list='S')
                            
                            # 瑞利波窗口
                            Rayleigh_begin = distance_in_km / Rayleigh_velocity[1]
                            Rayleigh_end = distance_in_km / Rayleigh_velocity[0]
                            
                            # 读取仅头部
                            sac = SACTrace.read(filename, headonly=True)
                            
                            # 写入SAC头部信息
                            sac.knetwk = network
                            sac.kstnm = station
                            sac.khole = location
                            sac.kcmpnm = channel
                            
                            sac.stla = sta_lat
                            sac.stlo = sta_lon
                            sac.stel = sta_ele
                            
                            sac.evla = evt_lat
                            sac.evlo = evt_lon
                            sac.evdp = evt_dep
                            sac.mag = magnitude
                            
                            sac.gcarc = distance_in_degree
                            sac.dist = distance_in_km
                            sac.az = az
                            sac.baz = baz
                            
                            sac.o = otime
                            sac.iztype = 'io'
                            sac.ko = 'O'
                            
                            # P波和S波到时
                            if P_arrivals:
                                sac.a = otime + P_arrivals[0].time
                                sac.t0 = otime + S_arrivals[0].time
                                sac.ka = 'P'
                                sac.kt0 = 'S'
                            
                            # 瑞利波窗口
                            sac.t1 = otime + Rayleigh_begin
                            sac.t2 = otime + Rayleigh_end
                            sac.kt1 = 'RayStart'
                            sac.kt2 = 'RayEnd'
                            
                            # 写入发震时刻到SAC头部
                            sac.kevnm = otime_str[0:4] + otime_str[5:7] + otime_str[8:10] + \
                                        otime_str[10:13] + otime_str[14:16] + otime_str[17:19]
                            # 发震时刻的小数部分
                            sac.kuser0 = otime_str[19:27]
                            # 震级类型
                            sac.kuser1 = magnitude_type
                            
                            # 只写入头部，文件必须存在
                            sac.write(filename, headonly=True)
                            
                            self.statusUpdate.emit(f'文件 {filename} 截取并保存成功')
                
                # 关闭进度条
                progress_adapter.close()
                self.statusUpdate.emit("事件数据截取完成！")
                
        except Exception as e:
            self.statusUpdate.emit(f"发生错误: {str(e)}")
        finally:
            self.finished.emit()


class SeismicDataDownloader(QMainWindow):
    def __init__(self):
        super().__init__()
        self.initUI()
        
    def initUI(self):
        self.setWindowTitle('地震数据下载工具')
        self.setGeometry(100, 100, 900, 700)
        
        # 创建工具栏
        self.create_toolbar()
        
        # 创建状态栏
        self.statusBar().showMessage('就绪')
        
        # 创建中央部件
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        # 主布局
        main_layout = QVBoxLayout(central_widget)
        
        # 创建选项卡部件
        self.tabs = QTabWidget()
        main_layout.addWidget(self.tabs)
        
        # 创建各个选项卡
        self.create_home_tab()  # 添加主页选项卡
        self.create_metadata_tab()
        self.create_daydata_tab()
        self.create_catalog_tab()
        self.create_eventdata_tab()
        self.create_cutevent_tab()
        
        # 状态和进度部分
        status_box = QGroupBox("状态和进度")
        status_layout = QVBoxLayout()
        
        self.statusText = QTextEdit()
        self.statusText.setReadOnly(True)
        self.progressBar = QProgressBar()
        
        status_layout.addWidget(self.statusText)
        status_layout.addWidget(self.progressBar)
        status_box.setLayout(status_layout)
        
        main_layout.addWidget(status_box)
        
    def create_toolbar(self):
        """创建工具栏"""
        self.toolbar = QToolBar("主工具栏")
        self.addToolBar(self.toolbar)
        
        # 打开数据目录操作
        open_data_action = QAction("打开数据目录", self)
        open_data_action.setStatusTip("打开存储数据的目录")
        open_data_action.triggered.connect(self.open_data_folder)
        self.toolbar.addAction(open_data_action)
        
        # 打开目录目录操作
        open_catalog_action = QAction("打开目录目录", self)
        open_catalog_action.setStatusTip("打开存储地震目录的目录")
        open_catalog_action.triggered.connect(self.open_catalog_folder)
        self.toolbar.addAction(open_catalog_action)
        
        self.toolbar.addSeparator()
        
        # 检查IRIS连接
        check_iris_action = QAction("检查IRIS连接", self)
        check_iris_action.setStatusTip("检查与IRIS数据中心的连接")
        check_iris_action.triggered.connect(self.check_iris_connection)
        self.toolbar.addAction(check_iris_action)
        
        self.toolbar.addSeparator()
        
        # 帮助操作
        help_action = QAction("帮助", self)
        help_action.setStatusTip("显示帮助信息")
        help_action.triggered.connect(lambda: self.tabs.setCurrentIndex(0))  # 切换到主页选项卡
        self.toolbar.addAction(help_action)
        
        # 关于操作
        about_action = QAction("关于", self)
        about_action.setStatusTip("关于此程序")
        about_action.triggered.connect(self.show_about)
        self.toolbar.addAction(about_action)
        
    def open_data_folder(self):
        """打开数据目录"""
        data_path = os.path.join(os.getcwd(), "DATA")
        if not os.path.exists(data_path):
            os.makedirs(data_path)
        QFileDialog.getOpenFileName(self, "打开数据目录", data_path)
        # 在Windows上尝试直接打开资源管理器
        os.startfile(data_path) if os.name == 'nt' else os.system(f'xdg-open "{data_path}"')
        
    def open_catalog_folder(self):
        """打开目录目录"""
        catalog_path = os.path.join(os.getcwd(), "catalog")
        if not os.path.exists(catalog_path):
            os.makedirs(catalog_path)
        QFileDialog.getOpenFileName(self, "打开目录目录", catalog_path)
        # 在Windows上尝试直接打开资源管理器
        os.startfile(catalog_path) if os.name == 'nt' else os.system(f'xdg-open "{catalog_path}"')
        
    def check_iris_connection(self):
        """检查与IRIS数据中心的连接"""
        try:
            client = Client('IRIS')
            QMessageBox.information(self, "连接状态", "与IRIS数据中心的连接正常。")
            self.update_status("IRIS连接检查成功！")
        except Exception as e:
            QMessageBox.critical(self, "连接状态", f"无法连接到IRIS数据中心：{str(e)}")
            self.update_status(f"IRIS连接检查失败：{str(e)}")
            
    def show_about(self):
        """显示关于对话框"""
        QMessageBox.about(self, "关于", 
                          "地震数据下载工具\n\n"
                          "版本：1.0\n"
                          "作者：吴悦初 (12131066@mail.sustech.edu.cn)\n\n"
                          "基于PyQt5和ObsPy开发，用于自动下载和处理地震数据。")
        
    def create_home_tab(self):
        """创建主页选项卡，提供使用指南和帮助信息"""
        tab = QWidget()
        layout = QVBoxLayout()
        
        # 欢迎标题
        welcome_label = QLabel("欢迎使用地震数据下载工具")
        welcome_label.setStyleSheet("font-size: 18pt; font-weight: bold;")
        welcome_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(welcome_label)
        
        # 使用指南
        guide_group = QGroupBox("使用指南")
        guide_layout = QVBoxLayout()
        
        guide_text = QTextEdit()
        guide_text.setReadOnly(True)
        guide_text.setHtml("""
        <p style='font-size: 11pt;'>本工具用于从IRIS数据中心下载和处理地震数据。请按照以下步骤操作：</p>
        <ol style='font-size: 11pt;'>
            <li><b>获取台站元数据</b>：首先打开"获取台站元数据"选项卡，设置网络代码、时间范围和通道查询参数，然后点击"获取台站元数据"按钮。</li>
            <li><b>获取地震目录</b>：接着打开"获取地震目录"选项卡，设置参数后点击"获取地震目录"按钮。</li>
            <li><b>下载数据</b>：您可以选择以下两种方式之一下载数据：
                <ul>
                    <li>打开"下载日数据"选项卡，下载连续的日数据。</li>
                    <li>打开"下载事件数据"选项卡，直接下载地震事件数据。</li>
                </ul>
            </li>
            <li><b>截取事件数据</b>：如果您已经下载了日数据，可以打开"截取事件数据"选项卡，从日数据中提取事件数据。</li>
        </ol>
        <p style='font-size: 11pt;'>在每个选项卡中，设置好参数后点击相应的按钮开始执行任务。下载和处理过程的进度和状态将显示在窗口底部。</p>
        """)
        guide_layout.addWidget(guide_text)
        guide_group.setLayout(guide_layout)
        layout.addWidget(guide_group)
        
        # 关于部分
        about_group = QGroupBox("关于")
        about_layout = QVBoxLayout()
        
        about_text = QTextEdit()
        about_text.setReadOnly(True)
        about_text.setHtml("""
        <p style='font-size: 11pt;'>作者：吴悦初 (12131066@mail.sustech.edu.cn)</p>
        <p style='font-size: 11pt;'>本工具基于PyQt5和ObsPy开发，用于自动下载和处理地震数据。</p>
        <p style='font-size: 11pt;'>如有问题或建议，请联系作者。</p>
        """)
        about_layout.addWidget(about_text)
        about_group.setLayout(about_layout)
        layout.addWidget(about_group)
        
        # 开始按钮
        start_button = QPushButton("开始使用 >>")
        start_button.setStyleSheet("font-size: 12pt; padding: 10px;")
        start_button.clicked.connect(lambda: self.tabs.setCurrentIndex(1))  # 点击后切换到元数据选项卡
        layout.addWidget(start_button)
        
        tab.setLayout(layout)
        self.tabs.addTab(tab, "主页")
        
    def create_metadata_tab(self):
        tab = QWidget()
        layout = QVBoxLayout()
        
        # 网络和时间设置
        network_group = QGroupBox("台站元数据设置")
        network_layout = QGridLayout()
        
        # 网络代码
        network_layout.addWidget(QLabel("网络代码:"), 0, 0)
        self.network_code_input = QLineEdit("XO")
        network_layout.addWidget(self.network_code_input, 0, 1)
        
        # 起始时间
        network_layout.addWidget(QLabel("起始时间:"), 1, 0)
        self.start_time_input = QLineEdit("2018-01-01")
        network_layout.addWidget(self.start_time_input, 1, 1)
        
        # 结束时间
        network_layout.addWidget(QLabel("结束时间:"), 2, 0)
        self.end_time_input = QLineEdit("2019-12-31")
        network_layout.addWidget(self.end_time_input, 2, 1)
        
        # 通道查询
        network_layout.addWidget(QLabel("通道查询:"), 3, 0)
        self.channel_query_input = QLineEdit("?H?,??H")
        network_layout.addWidget(self.channel_query_input, 3, 1)
        
        network_group.setLayout(network_layout)
        layout.addWidget(network_group)
        
        # 运行按钮
        run_button = QPushButton("获取台站元数据")
        run_button.clicked.connect(self.run_get_metadata)
        layout.addWidget(run_button)
        
        # 填充剩余空间
        layout.addStretch()
        
        tab.setLayout(layout)
        self.tabs.addTab(tab, "获取台站元数据")
        
    def create_daydata_tab(self):
        tab = QWidget()
        layout = QVBoxLayout()
        
        # 日数据设置
        daydata_group = QGroupBox("日数据下载设置")
        daydata_layout = QGridLayout()
        
        # 输出目录
        daydata_layout.addWidget(QLabel("输出目录:"), 0, 0)
        self.daydata_output_dir = QLineEdit("DATA/sacdata_day")
        browse_button = QPushButton("浏览...")
        browse_button.clicked.connect(lambda: self.browse_directory(self.daydata_output_dir))
        daydata_layout.addWidget(self.daydata_output_dir, 0, 1)
        daydata_layout.addWidget(browse_button, 0, 2)
        
        # 网络代码
        daydata_layout.addWidget(QLabel("网络代码:"), 1, 0)
        self.daydata_network = QLineEdit("XO")
        daydata_layout.addWidget(self.daydata_network, 1, 1)
        
        # 位置代码
        daydata_layout.addWidget(QLabel("位置代码:"), 2, 0)
        self.daydata_location = QLineEdit("--")
        daydata_layout.addWidget(self.daydata_location, 2, 1)
        
        # 台站列表
        daydata_layout.addWidget(QLabel("台站列表(逗号分隔):"), 3, 0)
        self.daydata_stations = QLineEdit("WD58,WS75")
        daydata_layout.addWidget(self.daydata_stations, 3, 1)
        
        # 地震计输出单位
        daydata_layout.addWidget(QLabel("地震计输出单位:"), 4, 0)
        self.daydata_units_seis = QComboBox()
        self.daydata_units_seis.addItems(["DISP", "VEL", "ACC"])
        daydata_layout.addWidget(self.daydata_units_seis, 4, 1)
        
        # 压力计输出单位
        daydata_layout.addWidget(QLabel("压力计输出单位:"), 5, 0)
        self.daydata_units_pres = QComboBox()
        self.daydata_units_pres.addItems(["DEF"])
        daydata_layout.addWidget(self.daydata_units_pres, 5, 1)
        
        # 元数据文件
        daydata_layout.addWidget(QLabel("元数据文件:"), 6, 0)
        self.daydata_metadata = QLineEdit("XO_fdsn_metadata.txt")
        browse_button2 = QPushButton("浏览...")
        browse_button2.clicked.connect(lambda: self.browse_file(self.daydata_metadata))
        daydata_layout.addWidget(self.daydata_metadata, 6, 1)
        daydata_layout.addWidget(browse_button2, 6, 2)
        
        daydata_group.setLayout(daydata_layout)
        layout.addWidget(daydata_group)
        
        # 运行按钮
        run_button = QPushButton("下载日数据")
        run_button.clicked.connect(self.run_download_daydata)
        layout.addWidget(run_button)
        
        # 填充剩余空间
        layout.addStretch()
        
        tab.setLayout(layout)
        self.tabs.addTab(tab, "下载日数据")
        
    def create_catalog_tab(self):
        tab = QWidget()
        layout = QVBoxLayout()
        
        # 地震目录设置
        catalog_group = QGroupBox("地震目录设置")
        catalog_layout = QGridLayout()
        
        # 输出目录
        catalog_layout.addWidget(QLabel("输出目录:"), 0, 0)
        self.catalog_output_dir = QLineEdit("catalog")
        browse_button = QPushButton("浏览...")
        browse_button.clicked.connect(lambda: self.browse_directory(self.catalog_output_dir))
        catalog_layout.addWidget(self.catalog_output_dir, 0, 1)
        catalog_layout.addWidget(browse_button, 0, 2)
        
        # 网络代码
        catalog_layout.addWidget(QLabel("网络代码:"), 1, 0)
        self.catalog_network = QLineEdit("XO")
        catalog_layout.addWidget(self.catalog_network, 1, 1)
        
        # 本地最小震级
        catalog_layout.addWidget(QLabel("本地最小震级:"), 2, 0)
        self.catalog_minmag_local = QLineEdit("5.5")
        catalog_layout.addWidget(self.catalog_minmag_local, 2, 1)
        
        # 本地最大震级
        catalog_layout.addWidget(QLabel("本地最大震级:"), 3, 0)
        self.catalog_maxmag_local = QLineEdit("6")
        catalog_layout.addWidget(self.catalog_maxmag_local, 3, 1)
        
        # 本地最大半径（度）
        catalog_layout.addWidget(QLabel("本地最大半径(度):"), 4, 0)
        self.catalog_maxradius_local = QLineEdit("90")
        catalog_layout.addWidget(self.catalog_maxradius_local, 4, 1)
        
        # 全球最小震级
        catalog_layout.addWidget(QLabel("全球最小震级:"), 5, 0)
        self.catalog_minmag_global = QLineEdit("6")
        catalog_layout.addWidget(self.catalog_minmag_global, 5, 1)
        
        # 元数据文件
        catalog_layout.addWidget(QLabel("元数据文件:"), 6, 0)
        self.catalog_metadata = QLineEdit("XO_fdsn_metadata.txt")
        browse_button2 = QPushButton("浏览...")
        browse_button2.clicked.connect(lambda: self.browse_file(self.catalog_metadata))
        catalog_layout.addWidget(self.catalog_metadata, 6, 1)
        catalog_layout.addWidget(browse_button2, 6, 2)
        
        catalog_group.setLayout(catalog_layout)
        layout.addWidget(catalog_group)
        
        # 运行按钮
        run_button = QPushButton("获取地震目录")
        run_button.clicked.connect(self.run_get_catalog)
        layout.addWidget(run_button)
        
        # 填充剩余空间
        layout.addStretch()
        
        tab.setLayout(layout)
        self.tabs.addTab(tab, "获取地震目录")
        
    def create_eventdata_tab(self):
        tab = QWidget()
        layout = QVBoxLayout()
        
        # 事件数据设置
        event_group = QGroupBox("事件数据下载设置")
        event_layout = QGridLayout()
        
        # 输出目录
        event_layout.addWidget(QLabel("输出目录:"), 0, 0)
        self.event_output_dir = QLineEdit("DATA/sacdata_event")
        browse_button = QPushButton("浏览...")
        browse_button.clicked.connect(lambda: self.browse_directory(self.event_output_dir))
        event_layout.addWidget(self.event_output_dir, 0, 1)
        event_layout.addWidget(browse_button, 0, 2)
        
        # 网络代码
        event_layout.addWidget(QLabel("网络代码:"), 1, 0)
        self.event_network = QLineEdit("XO")
        event_layout.addWidget(self.event_network, 1, 1)
        
        # 位置代码
        event_layout.addWidget(QLabel("位置代码:"), 2, 0)
        self.event_location = QLineEdit("--")
        event_layout.addWidget(self.event_location, 2, 1)
        
        # 台站列表
        event_layout.addWidget(QLabel("台站列表(逗号分隔):"), 3, 0)
        self.event_stations = QLineEdit("WD58,WS75")
        event_layout.addWidget(self.event_stations, 3, 1)
        
        # 地震计输出单位
        event_layout.addWidget(QLabel("地震计输出单位:"), 4, 0)
        self.event_units_seis = QComboBox()
        self.event_units_seis.addItems(["DISP", "VEL", "ACC"])
        event_layout.addWidget(self.event_units_seis, 4, 1)
        
        # 压力计输出单位
        event_layout.addWidget(QLabel("压力计输出单位:"), 5, 0)
        self.event_units_pres = QComboBox()
        self.event_units_pres.addItems(["DEF"])
        event_layout.addWidget(self.event_units_pres, 5, 1)
        
        # 事件长度
        event_layout.addWidget(QLabel("事件长度(秒):"), 6, 0)
        self.event_length = QLineEdit("7200")
        event_layout.addWidget(self.event_length, 6, 1)
        
        # 起始时间偏移
        event_layout.addWidget(QLabel("起始时间偏移(秒):"), 7, 0)
        self.event_btime = QLineEdit("-600")
        event_layout.addWidget(self.event_btime, 7, 1)
        
        # 瑞利波速度范围
        event_layout.addWidget(QLabel("瑞利波速度范围(km/s):"), 8, 0)
        self.event_rayleigh_vel = QLineEdit("3,4.5")
        event_layout.addWidget(self.event_rayleigh_vel, 8, 1)
        
        # 元数据文件
        event_layout.addWidget(QLabel("元数据文件:"), 9, 0)
        self.event_metadata = QLineEdit("XO_fdsn_metadata.txt")
        browse_button2 = QPushButton("浏览...")
        browse_button2.clicked.connect(lambda: self.browse_file(self.event_metadata))
        event_layout.addWidget(self.event_metadata, 9, 1)
        event_layout.addWidget(browse_button2, 9, 2)
        
        event_group.setLayout(event_layout)
        layout.addWidget(event_group)
        
        # 运行按钮
        run_button = QPushButton("下载事件数据")
        run_button.clicked.connect(self.run_download_eventdata)
        layout.addWidget(run_button)
        
        # 填充剩余空间
        layout.addStretch()
        
        tab.setLayout(layout)
        self.tabs.addTab(tab, "下载事件数据")
        
    def create_cutevent_tab(self):
        tab = QWidget()
        layout = QVBoxLayout()
        
        # 截取事件设置
        cut_group = QGroupBox("从日数据中截取事件数据设置")
        cut_layout = QGridLayout()
        
        # 输入目录
        cut_layout.addWidget(QLabel("输入目录:"), 0, 0)
        self.cut_input_dir = QLineEdit("DATA/sacdata_day")
        browse_button = QPushButton("浏览...")
        browse_button.clicked.connect(lambda: self.browse_directory(self.cut_input_dir))
        cut_layout.addWidget(self.cut_input_dir, 0, 1)
        cut_layout.addWidget(browse_button, 0, 2)
        
        # 输出目录
        cut_layout.addWidget(QLabel("输出目录:"), 1, 0)
        self.cut_output_dir = QLineEdit("DATA/sacdata_event")
        browse_button2 = QPushButton("浏览...")
        browse_button2.clicked.connect(lambda: self.browse_directory(self.cut_output_dir))
        cut_layout.addWidget(self.cut_output_dir, 1, 1)
        cut_layout.addWidget(browse_button2, 1, 2)
        
        # 目录目录
        cut_layout.addWidget(QLabel("目录目录:"), 2, 0)
        self.cut_catalog_dir = QLineEdit("catalog")
        browse_button3 = QPushButton("浏览...")
        browse_button3.clicked.connect(lambda: self.browse_directory(self.cut_catalog_dir))
        cut_layout.addWidget(self.cut_catalog_dir, 2, 1)
        cut_layout.addWidget(browse_button3, 2, 2)
        
        # 网络代码
        cut_layout.addWidget(QLabel("网络代码:"), 3, 0)
        self.cut_network = QLineEdit("XO")
        cut_layout.addWidget(self.cut_network, 3, 1)
        
        # 位置代码
        cut_layout.addWidget(QLabel("位置代码:"), 4, 0)
        self.cut_location = QLineEdit("--")
        cut_layout.addWidget(self.cut_location, 4, 1)
        
        # 台站列表
        cut_layout.addWidget(QLabel("台站列表(逗号分隔):"), 5, 0)
        self.cut_stations = QLineEdit("WD58,WS75")
        cut_layout.addWidget(self.cut_stations, 5, 1)
        
        # 事件长度
        cut_layout.addWidget(QLabel("事件长度(秒):"), 6, 0)
        self.cut_event_length = QLineEdit("7200")
        cut_layout.addWidget(self.cut_event_length, 6, 1)
        
        # 起始时间偏移
        cut_layout.addWidget(QLabel("起始时间偏移(秒):"), 7, 0)
        self.cut_btime = QLineEdit("-600")
        cut_layout.addWidget(self.cut_btime, 7, 1)
        
        # 瑞利波速度范围
        cut_layout.addWidget(QLabel("瑞利波速度范围(km/s):"), 8, 0)
        self.cut_rayleigh_vel = QLineEdit("3,4.5")
        cut_layout.addWidget(self.cut_rayleigh_vel, 8, 1)
        
        # 元数据文件
        cut_layout.addWidget(QLabel("元数据文件:"), 9, 0)
        self.cut_metadata = QLineEdit("XO_fdsn_metadata.txt")
        browse_button4 = QPushButton("浏览...")
        browse_button4.clicked.connect(lambda: self.browse_file(self.cut_metadata))
        cut_layout.addWidget(self.cut_metadata, 9, 1)
        cut_layout.addWidget(browse_button4, 9, 2)
        
        cut_group.setLayout(cut_layout)
        layout.addWidget(cut_group)
        
        # 运行按钮
        run_button = QPushButton("截取事件数据")
        run_button.clicked.connect(self.run_cut_eventdata)
        layout.addWidget(run_button)
        
        # 填充剩余空间
        layout.addStretch()
        
        tab.setLayout(layout)
        self.tabs.addTab(tab, "截取事件数据")
        
    def browse_directory(self, line_edit):
        directory = QFileDialog.getExistingDirectory(self, "选择目录")
        if directory:
            line_edit.setText(directory)
            
    def browse_file(self, line_edit):
        file, _ = QFileDialog.getOpenFileName(self, "选择文件", "", "文本文件 (*.txt);;所有文件 (*.*)")
        if file:
            line_edit.setText(file)
    
    def update_status(self, message):
        self.statusText.append(message)
        # 滚动到底部
        self.statusText.verticalScrollBar().setValue(self.statusText.verticalScrollBar().maximum())
        # 更新状态栏
        self.statusBar().showMessage(message)
        
    def update_progress(self, value):
        self.progressBar.setValue(value)
        
    def run_get_metadata(self):
        # 显示确认对话框
        reply = QMessageBox.question(self, '确认操作', 
                                      '即将开始获取台站元数据，这个过程可能需要几分钟时间。\n\n确定要继续吗？',
                                      QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
        
        if reply == QMessageBox.No:
            return
            
        # 从UI获取参数
        from obspy import UTCDateTime
        params = {
            "network_code": self.network_code_input.text(),
            "network_starttime": UTCDateTime(self.start_time_input.text()),
            "network_endtime": UTCDateTime(self.end_time_input.text()),
            "channel_query": self.channel_query_input.text()
        }
        
        # 创建工作线程
        self.worker = WorkerThread("metadata", params)
        self.worker.statusUpdate.connect(self.update_status)
        self.worker.progressUpdate.connect(self.update_progress)
        self.worker.finished.connect(self.on_worker_finished)
        
        # 开始执行
        self.update_status("开始获取台站元数据...")
        self.worker.start()
        
    def run_download_daydata(self):
        # 显示确认对话框
        reply = QMessageBox.question(self, '确认操作', 
                                      '即将开始下载日数据，这个过程可能需要较长时间（取决于数据量大小）。\n\n确定要继续吗？',
                                      QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
        
        if reply == QMessageBox.No:
            return
            
        # 从UI获取参数
        params = {
            "output_dir": self.daydata_output_dir.text(),
            "network": self.daydata_network.text(),
            "location": self.daydata_location.text(),
            "stations_download": self.daydata_stations.text().split(','),
            "output_units_seis": self.daydata_units_seis.currentText(),
            "output_units_pres": self.daydata_units_pres.currentText(),
            "metadatafile": self.daydata_metadata.text()
        }
        
        # 创建工作线程
        self.worker = WorkerThread("daydata", params)
        self.worker.statusUpdate.connect(self.update_status)
        self.worker.progressUpdate.connect(self.update_progress)
        self.worker.finished.connect(self.on_worker_finished)
        
        # 开始执行
        self.update_status("开始下载日数据...")
        self.worker.start()
        
    def run_get_catalog(self):
        # 显示确认对话框
        reply = QMessageBox.question(self, '确认操作', 
                                      '即将开始获取地震目录，这个过程可能需要几分钟时间。\n\n确定要继续吗？',
                                      QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
        
        if reply == QMessageBox.No:
            return
            
        # 从UI获取参数
        params = {
            "output_dir": self.catalog_output_dir.text(),
            "network": self.catalog_network.text(),
            "minmag_local": float(self.catalog_minmag_local.text()),
            "maxmag_local": float(self.catalog_maxmag_local.text()),
            "maxradius_local": float(self.catalog_maxradius_local.text()),
            "minmag_global": float(self.catalog_minmag_global.text()),
            "metadatafile": self.catalog_metadata.text()
        }
        
        # 创建工作线程
        self.worker = WorkerThread("catalog", params)
        self.worker.statusUpdate.connect(self.update_status)
        self.worker.progressUpdate.connect(self.update_progress)
        self.worker.finished.connect(self.on_worker_finished)
        
        # 开始执行
        self.update_status("开始获取地震目录...")
        self.worker.start()
        
    def run_download_eventdata(self):
        # 显示确认对话框
        reply = QMessageBox.question(self, '确认操作', 
                                      '即将开始下载事件数据，这个过程可能需要较长时间（取决于事件数量）。\n\n确定要继续吗？',
                                      QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
        
        if reply == QMessageBox.No:
            return
            
        # 从UI获取参数
        rayleigh_vel = [float(x) for x in self.event_rayleigh_vel.text().split(',')]
        params = {
            "output_dir": self.event_output_dir.text(),
            "network": self.event_network.text(),
            "location": self.event_location.text(),
            "stations_download": self.event_stations.text().split(','),
            "output_units_seis": self.event_units_seis.currentText(),
            "output_units_pres": self.event_units_pres.currentText(),
            "event_length": float(self.event_length.text()),
            "btime": float(self.event_btime.text()),
            "Rayleigh_velocity": rayleigh_vel,
            "metadatafile": self.event_metadata.text()
        }
        
        # 创建工作线程
        self.worker = WorkerThread("eventdata", params)
        self.worker.statusUpdate.connect(self.update_status)
        self.worker.progressUpdate.connect(self.update_progress)
        self.worker.finished.connect(self.on_worker_finished)
        
        # 开始执行
        self.update_status("开始下载事件数据...")
        self.worker.start()
        
    def run_cut_eventdata(self):
        # 显示确认对话框
        reply = QMessageBox.question(self, '确认操作', 
                                      '即将开始从日数据中截取事件数据，这个过程可能需要一些时间。\n\n确定要继续吗？',
                                      QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
        
        if reply == QMessageBox.No:
            return
            
        # 从UI获取参数
        rayleigh_vel = [float(x) for x in self.cut_rayleigh_vel.text().split(',')]
        params = {
            "input_dir": self.cut_input_dir.text(),
            "output_dir": self.cut_output_dir.text(),
            "catlog_dir": self.cut_catalog_dir.text(),
            "network": self.cut_network.text(),
            "location": self.cut_location.text(),
            "stations_download": self.cut_stations.text().split(','),
            "event_length": float(self.cut_event_length.text()),
            "btime": float(self.cut_btime.text()),
            "Rayleigh_velocity": rayleigh_vel,
            "metadatafile": self.cut_metadata.text()
        }
        
        # 创建工作线程
        self.worker = WorkerThread("cutevents", params)
        self.worker.statusUpdate.connect(self.update_status)
        self.worker.progressUpdate.connect(self.update_progress)
        self.worker.finished.connect(self.on_worker_finished)
        
        # 开始执行
        self.update_status("开始从日数据中截取事件数据...")
        self.worker.start()
        
    def on_worker_finished(self):
        self.update_status("任务完成！")
        # 更新状态栏
        self.statusBar().showMessage("任务完成 - 就绪")
        QMessageBox.information(self, "完成", "操作已完成！")


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = SeismicDataDownloader()
    # 显示主界面
    ex.show()
    # 在这里显示一条欢迎消息，告诉用户如何开始
    ex.update_status("欢迎使用地震数据下载工具！请选择相应的选项卡，设置参数后点击对应的按钮开始下载。")
    sys.exit(app.exec_()) 