# Mendelian-Randomization
In this repository,I will start a project to show how to do mendelian randomization

# 血液中物质GWAS数据的获取

首先在GWAScatalog网站上检索文章Genetics of 35 blood and urine biomarkers in the UK Biobank的GWAS数据

![image](https://user-images.githubusercontent.com/105405121/174446161-7098d670-408b-480d-8544-c4ea505b9efd.png)

点击检索结果，进入检索结果页面，划到最下方，如下图所示：

![image](https://user-images.githubusercontent.com/105405121/174446695-2004020f-ea10-4672-b7c5-0f78e17efe07.png)

上图中左侧的Study accession一列即为35种物质的GWAS数据列表，点击其中任何一个，即可看到如下页面：

![image](https://user-images.githubusercontent.com/105405121/174446851-44324de3-b453-438c-9f69-78a2ab7198d5.png)

点击上图中的FTP Download链接，结果如下：

![image](https://user-images.githubusercontent.com/105405121/174446887-44a91bbd-4e7e-4ea9-b152-366f902b2cba.png)

上图中以.tsv结尾的文件即GWAS数据文件，也就是以制表符分隔的文本文件

# covid19 GWAS数据的获取

进入https://www.covid19hg.org/ ，点击左侧工具栏中Results-Downloads-Release6：

![image](https://user-images.githubusercontent.com/105405121/174447117-b99241a7-7c1d-4f6d-bb67-f25c3a272f02.png)

向下翻动即可看到Download相关链接，选择38版本的数据进行下载：

![image](https://user-images.githubusercontent.com/105405121/174447145-bf699c33-c802-4a3a-9011-0289938a41bf.png)

# 用TwoSampleMR包进行孟德尔随机化

library(TwoSampleMR)//载入R包
setwd("D:/GWAS文件/35/性激素结合球蛋白")//将工作路径设置为GWAS数据文件所在的目录下
danbai <- read.table("GCST90019518_buildGRCh37.tsv",header = T)//GCST90019518_buildGRCh37.tsv为血液中物质的GWAS数据(1/35)
exp_dat <- format_data(danbai,type = "exposure",snp_col = "variant_id",beta_col = "beta",se_col = "standard_error",effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p_value")
danbai_exp_dat <- clump_data(exp_dat,clump_kb = 1000000,clump_r2 = 0.00001)//去除连锁不平衡
library(data.table)
covid19 <- fread("COVID19_HGI_C2_ALL_leave_23andme_20210607.txt",header = T)
covid19_out_dat <- format_data(covid19,type = "outcome",snps = danbai_exp_dat$SNP,header = TRUE,snp_col = "rsid",beta_col = "all_inv_var_meta_beta",se_col = "all_inv_var_meta_sebeta",effect_allele_col = "ALT",other_allele_col = "REF",pval_col = "all_inv_var_meta_p",ncase_col = "all_inv_var_meta_cases",ncontrol_col = "all_inv_var_meta_controls",chr_col = "#CHR",pos_col = "POS")//COVID19 GWAS数据每一列的含义可以在之前的covidhg.org网站上找到对应的说明
mydata <- harmonise_data(exposure_dat = danbai_exp_dat,outcome_dat = covid19_out_dat,action = 2)//将IV的效应等位基因对齐
res <- mr(mydata)//孟德尔随机化
het <- mr_heterogeneity(mydata)//异质性检验
pleio <- mr_pleiotropy_test(mydata)//水平多样性检验
single <- mr_leaveoneout(mydata)
mr_leaveoneout_plot(single)//逐一剔除检验
p1 <- mr_scatter_plot(res,mydata)//散点图
res_single <- mr_singlesnp(mydata)
p2 <- mr_forest_plot(res_single)//森林图
p3 <- mr_funnel_plot(res_single)//漏斗图
p1[[1]]//查看散点图
p2[[1]]//查看森林图
p3[[1]]//查看漏斗图
