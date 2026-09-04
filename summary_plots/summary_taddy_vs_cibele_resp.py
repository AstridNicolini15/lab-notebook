#%%
folders_cibele = ["SST-cells_WT_Adult_V1",
    "SST-cells_cond-GluN1-KO_Adult_V1"]
summary_path_cibele = ["/home/user/DATA/Astrid/Cibele_data/summary"] * len(folders_cibele)
colors_cibele = [[(to_rgb('#12522eff'), 0.8),'lightgrey'], [(to_rgb('#5d1490ff'), 0.8),'lightgrey']]


folders_taddy = ["Wild-Type",
    "GluN1-KO"
    ]
summary_path_taddy = ["/home/user/DATA/Astrid/Taddy_data/OneDrive_1_9-3-2026/summary"] * len(folders_taddy)
colors_taddy = [[pt.tab10(1),'lightgrey'], [pt.tab10(4),'lightgrey']]


#%%

mosaic = """
    AAAAAAhBBBBBB
    AAAAAAhBBBBBB
    AAAAAAhBBBBBB
    AAAAAAhBBBBBB
    AAAAAAhBBBBBB
    ooooooooooooo
    EEFFjIIJJkCCC
    EEFFjIIJJkCCC
    KKLLjMMNNkCCC
    KKLLjMMNNkCCC
    """

fig = plt.figure(layout="constrained", figsize = (16,12))
ax_dict = fig.subplot_mosaic(mosaic)

ax_dict["h"].axis("off")
ax_dict["j"].axis("off")
ax_dict["o"].axis("off")
ax_dict["k"].axis("off")


for k, ax_letters, folders, summary_path, colors in zip([0,1],
                                                        [["A", "E", "F","K", "L"], ["B", "I", "J", "M", "N"]], 
                                                        [folders_cibele, folders_taddy], 
                                                        [summary_path_cibele, summary_path_taddy], 
                                                        [colors_cibele, colors_taddy]): 


    for i, folder in enumerate(folders) : 
        fig.text(x= 0+0.45*k, y= -0.02-0.04*i, s = '%s' % folder, color = colors[i][0], fontsize = 15)
    
    plot_tuning_responses_many_pop(folders, colors, ax_dict[ax_letters[0]], special_dict =  {'title' : '', 'ylims' : (-0.05,1.2)}, summary_path = summary_path)

    axes = [ax_dict[ax_letters[1]], ax_dict[ax_letters[3]], ax_dict[ax_letters[2]], ax_dict[ax_letters[4]]]
    plot_responsiveness_pie(folders, colors, axes, summary_path = summary_path)


folders = folders_cibele + folders_taddy
summary_path = summary_path_cibele + summary_path_taddy
colors = colors_cibele + colors_taddy

colors_without_grey = []
for color in colors : 
    colors_without_grey.append(color[0])
plot_mean_F_val(folders, 'correctedFluo0', colors_without_grey, ax_dict["C"], summary_path = summary_path)
