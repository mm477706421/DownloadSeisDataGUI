"""
进度条适配器 - 将tqdm进度条的更新转换为PyQt信号

Yuechu Wu
12131066@mail.sustech.edu.cn
2024-05-06
"""

class TqdmProgressAdapter:
    """
    tqdm进度条适配器，用于将tqdm进度条的更新转换为PyQt信号
    """
    def __init__(self, signal_progress=None, signal_status=None, total=100):
        self.signal_progress = signal_progress
        self.signal_status = signal_status
        self.total = total
        self.n = 0
        self.last_print_n = 0
        
    def update(self, n=1):
        """更新进度"""
        self.n += n
        percentage = int(100 * self.n / self.total) if self.total else 0
        
        # 发送进度信号
        if self.signal_progress:
            self.signal_progress.emit(percentage)
            
        # 当更新足够时发送状态信号
        if self.signal_status and (self.n - self.last_print_n) > (self.total / 100):
            self.signal_status.emit(f"进度: {self.n}/{self.total} ({percentage}%)")
            self.last_print_n = self.n
            
    def close(self):
        """关闭进度条"""
        if self.signal_progress:
            self.signal_progress.emit(100)
            
        if self.signal_status:
            self.signal_status.emit(f"进度: {self.total}/{self.total} (100%)") 