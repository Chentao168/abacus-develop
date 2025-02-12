from ase.io.cube import read_cube_data
import matplotlib.pyplot as plt
import numpy as np

# 读取CUBE文件
cube_data, atoms = read_cube_data('OUT.autotest/SPIN1_LDOS_1eV.cube')
# cube_data, atoms = read_cube_data('OUT.autotest/SPIN1_CHG.cube')

print(cube_data.shape)
print(atoms)

m, n = 3, 5
z = 8

z_index = z / atoms.cell[2, 2]
print(z_index)
z_index = int(z / atoms.cell[2, 2] * cube_data.shape[2])
print(z_index)

# 计算cube_data在z上的值
cube_data_z = cube_data[:, :, z_index]
cube_data_z = np.tile(cube_data_z, (m, n))

# cube_data_z = np.sum(cube_data[:, :, z_index:], axis=2)
# cube_data_z = np.tile(cube_data_z, (m, n))

# cube_data_z = np.sum(cube_data[:, :, : z_index + 1], axis=2)
# cube_data_z = np.tile(cube_data_z, (m, n))

# 生成下x, y的坐标
x = np.linspace(0, 1, cube_data.shape[0] * m)
y = np.linspace(0, 1, cube_data.shape[1] * n)
X, Y = np.meshgrid(x, y, indexing='ij')
# print(x)
# print(y)

r_x = X * atoms.cell[0, 0] * m + Y * atoms.cell[1, 0] * n
r_y = X * atoms.cell[0, 1] * m + Y * atoms.cell[1, 1] * n
# print(r_x)
# print(r_y)

# 画图

plt.figure(figsize=(8, 6))
plt.pcolormesh(r_x, r_y, cube_data_z, shading='auto', cmap='viridis')
plt.xlim(0, 16)
plt.ylim(-1, 13)
plt.colorbar(label='LDOS')
plt.xlabel('X (Å)')
plt.ylabel('Y (Å)')
plt.title('LDOS')
plt.savefig('ldos.png')
