
repeat.anova = function(data,n.g,n.t,N){

	if (!requireNamespace("tidyverse", quietly = TRUE))
	    install.packages("tidyverse")
	if (!requireNamespace("ggpubr", quietly = TRUE))
		install.packages("ggpubr")
	if (!requireNamespace("rstatix", quietly = TRUE))
	    install.packages("rstatix")
	if (!requireNamespace("datarium", quietly = TRUE))
		install.packages("datarium")
	if (!requireNamespace("do", quietly = TRUE))
		install.packages("do")
	library("tidyverse")
	library("ggpubr")
	library("rstatix")
	library("datarium")
	library("do")
	selfesteem2 = data.frame(id = seq(1:N),data)
	colnames(selfesteem2) = c("id","group",paste("t",seq(1:n.t),sep=""))
	selfesteem2 %>% sample_n_by(group, size = 1)

	selfesteem2 <- selfesteem2 %>%
		gather(key = "time", value = "score", paste("t",seq(1:n.t),sep="")) %>%
		convert_as_factor(id, time)

	selfesteem2 %>%
		group_by(group, time) %>%
		get_summary_stats(score, type = "mean_sd") ->c1
	  
	grouP = paste("G",seq(1:n.g),sep="")
	c.1 = NULL
	for(j in 1:n.g){
	c.tmp = c1[which(c1$group==grouP[j]),]
	c.tmp = paste(round(c.tmp$mean,2),"+-",round(c.tmp$sd,2))
	c.1 = rbind(c.1,c.tmp)
	}

	res.aov <- anova_test(
		data = selfesteem2, dv = score, wid = id,
		between = group, within = time
	)
	get_anova_table(res.aov) -> total.effect


	# Effect of treatment at each time point
	one.way <- selfesteem2 %>%
			group_by(time) %>%
			anova_test(dv = score, wid = id, between = group) %>%
			get_anova_table() %>%
			adjust_pvalue(method = "bonferroni")
	one.way -> t.effect
	c2 = matrix(NA,2,n.t+2)
	c2[,1:n.t] = rbind(t.effect$F,t.effect$p.adj)

	# Pairwise comparisons between group levels
	pwc.t <- selfesteem2 %>%
		group_by(time) %>%
		pairwise_t_test(score ~ group, p.adjust.method = "bonferroni")
	pwc.t

	# Effect of time at each level of exercises group
	one.way2 <- selfesteem2 %>%
			group_by(group) %>%
			anova_test(dv = score, wid = id, within = time) %>%
			get_anova_table() %>%
			adjust_pvalue(method = "bonferroni")
	one.way2 -> g.effect
	c3 = cbind(g.effect$F,g.effect$p.adj)

	# Pairwise comparisons between time points at each group levels
	# Paired t-test is used because we have repeated measures by time
	pwc.g <- selfesteem2 %>%
		group_by(group) %>%
		pairwise_t_test(
		score ~ time, paired = TRUE, 
		p.adjust.method = "bonferroni"
		) %>%
		select(-df, -statistic, -p) # Remove details
	pwc.g

	# comparisons for time variable
	selfesteem2 %>%
		pairwise_t_test(
		score ~ time, paired = TRUE, 
		p.adjust.method = "bonferroni"
		)
	# comparisons for group variable
	selfesteem2 %>%
		pairwise_t_test(
		score ~ group, 
		p.adjust.method = "bonferroni"
		)
	res = NULL
	res = cbind(c.1,c3)
	res = rbind(res,c2)
	colnames(res)  =c(paste("t",seq(1:n.t),sep=""),"F","P")
	rownames(res)  =c(grouP,"F","P")
	#res = data.frame(res)
	#res = apply(res,2,as.character)

	ref.a = data.frame(group = c("G1","G2","G3","G4","G5","G6","G7","G8"),sig = c("#","*","&","$","^","α","β","γ"),ind = seq(1:8))
	ref.b = data.frame(time = c("t1","t2","t3","t4","t5","t6","t7","t8"),sig = c("a","b","c","d","e","f","g","h"),ind = seq(1:8))

	n.1 = dim(pwc.t)
	for(i in 1:n.1){
		try({
		tmp = pwc.t[i,]
		if(tmp$p.adj<0.05){
			t.tmp = tmp$time
			g.tmp.1 = tmp$group1
			g.tmp.2 = tmp$group2
			sig.tmp = ref.a[which(ref.a$group==g.tmp.1),2]
			index.g = ref.a[which(ref.a$group==g.tmp.2),3]
			index.t = ref.b[which(ref.b$time==t.tmp),3]	
			res[index.g,index.t] <- as.character(paste(res[index.g,index.t],sig.tmp,sep=""))
	  }
	  })
	}

	n.2 = dim(pwc.g)
	for(i in 1:n.2){
		try({
		tmp = pwc.g[i,]
		if(tmp$p.adj<0.05){
			g.tmp = tmp$group
			t.tmp.1 = tmp$group1
			t.tmp.2 = tmp$group2
			sig.tmp = ref.b[which(ref.b$time==t.tmp.1),2]
			index.t = ref.b[which(ref.b$time==t.tmp.2),3]
			index.g = ref.a[which(ref.a$group==g.tmp),3]	
			res[index.g,index.t] <- as.character(paste(res[index.g,index.t],sig.tmp,sep=""))
		}
		})
	}

	res2 = paste("F:group = ",total.effect$F[1],", p:group = ",total.effect$p[1],
			", F:time = ",total.effect$F[2],", p:time = ",total.effect$p[2], 
			", F:interaction = ",total.effect$F[2],", p:interaction = ",total.effect$p[2],sep = "")
	ref.note = NULL
	for(z in 1:(n.t-1)){ref.tmp = paste(ref.b[z,2],":",ref.b[z,1],sep = "")
			ref.note = paste(ref.note,ref.tmp,sep="; ")}
	for(z in 1:(n.g-1)){ref.tmp = paste(ref.a[z,2],":",ref.a[z,1],sep = "")
			ref.note = paste(ref.note,ref.tmp,sep="; ")}
	return(list(res = res,effect = res2,note = ref.note))
	}		
