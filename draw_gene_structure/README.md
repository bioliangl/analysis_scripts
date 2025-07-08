*This project is a modified version of [wu116/Wu-self-analysis/pipeline/plot_gene_structure], created by [wu116] (GitHub: [wu116]). Original repository: [https://github.com/wu116/Wu-self-analysis]*

根据gff注释文件和interpro结构域注释，绘制转录本结构图

1、准备某个基因的注释文件

``/example/AthFCA.gff`` 拟南芥FCA基因，包含四个转录本

2、interpro (https://www.ebi.ac.uk/interpro/)注释结构域

interpro根据蛋白序列``/example/AthFCA.pep`` 得到注释结果，下载**GFF output** ``analysis/AthFCA.anno`` 

3、得到最后画图所需格式

```sh
python ./scripts/change_format.py -g /example/AthFCA.gff -a analysis/AthFCA.anno -d Pfam -s true -o analysis/input.gff
```

> **-d** 根据不同的数据库选择可视化的结构域，默认为："Pfam" 可选："Pfam", "CDD", "SUPERFAMILY"
>
> **-s** 是否标准化坐标显示(mRNA坐标从1开始)， 默认为："ture"， 显示原有位置则修改为"false"

4、R语言绘制结构图

```sh
Rscript draw_gene_structure.R -i analysis/input.gff -o /analysis/output -f png
```

> **-f** 输出的图片格式，默认为png，可选png、jpg、svg、pdf、pptx等，如果选择pptx需要保证已安装**R包eoffice** 

结果``./analysis/output.png``

![output](.\analysis\output.png)

