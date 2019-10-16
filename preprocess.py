import pandas as pd

# collecting trajectory data
col_names = ["x1","x2","x3","x4","x5","x6","x7","y1","y2","y3","y4","y5","y6","y7",
            "vx1","vx2","vx3","vx4","vx5","vx6","vx7","vy1","vy2","vy3","vy4","vy5","vy6","vy7",
            "ax1","ax2","ax3","ax4","ax5","ax6","ax7","ay1","ay2","ay3","ay4","ay5","ay6","ay7","J"]
df = pd.read_csv("/Users/admin/Desktop/V-formation/AMPC/traj.csv", names = col_names)

# relative position
for i in range(1,8):
    df['x' '%d' %(i)] = df['x' '%d' %(i)] - df['x1']
    df['y' '%d' %(i)] = df['y' '%d' %(i)] - df['y1']

# verifying data
print(df.head())

# saving dataframe
df.to_csv("/Users/admin/Desktop/V-formation/AMPC/newtraj.csv", index=False)
