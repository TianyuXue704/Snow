运行run_f_read_ATL09_strong.m，通过修改路径，可以批量读取指定路径下的所有ATL03*.h5文件中的强光束信息。
	其中调用了函数f_read_ATL09_strong.m
		其中调用了函数read_ATL03_strong1.m，read_ATL03_strong2.m，read_ATL03_strong3.m
		可以实现强光束1/2/3的光束信息读取
		由于我这里只对强光束3进行了处理，所以对函数read_ATL03_strong1.m，read_ATL03_strong2.m部分进行了注释。

运行run_f_read_ATL03_weak.m，同理，用于批量读取指定路径下的所有ATL03*.h5文件中的弱光束信息。
	有需要可以自己使用