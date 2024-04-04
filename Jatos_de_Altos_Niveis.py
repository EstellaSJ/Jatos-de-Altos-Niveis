# Instalando bibliotecas
!pip install cartopy
!pip install ncBuilder
!pip install matplotlib

# Importando as bibliotecas
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as mcolors
from datetime import datetime
from matplotlib.ticker import MultipleLocator
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER

# Carregando o arquivo com os dados de vento, pressão e altura geopotencial
vento250 = xr.open_dataset('vento.nc')
pressao = xr.open_dataset('pressao.nc')
geo = xr.open_dataset('geopotencial.nc')

vento_u = vento250['u']
vento_v = vento250['v']
latitude = vento250['latitude']
longitude = vento250['longitude']
pressao_mmsl_data = pressao['msl']
geopot = geo.variables['z'][:]

# Definindo as coordenadas
latitude_min = -60
latitude_max = 6
longitude_min = -75
longitude_max = -10

# Selecionando o passo de tempo
tempos = list(range(0, 160, 1))

for i in tempos:
    vento_u_t0 = vento_u.isel(time=i)
    vento_v_t0 = vento_v.isel(time=i)

    # Calculando a magnitude do vetor de vento
    magnitude = np.sqrt(vento_u_t0**2 + vento_v_t0**2)

    # Selecione um tempo específico (você pode ajustar conforme necessário)
    tempo_selecionado = 0
    pressao_mmsl_t0 = pressao_mmsl_data.isel(time=tempo_selecionado)

    # Criando a figura e adicionando as linhas do mapa
    fig, ax = plt.subplots(figsize=(15, 10), subplot_kw={'projection': ccrs.PlateCarree()})

    # Adicionando as fronteiras dos estados brasileiros
    ax.add_feature(cfeature.BORDERS.with_scale('10m'), linestyle='-', edgecolor='dimgray')
    ax.add_feature(cfeature.STATES.with_scale('10m'), linestyle='-', edgecolor='dimgray')

    # Subamostragem para reduzir a quantidade de vetores
    subsample_factor = 10  # Ajuste conforme necessário
    subsample_factor_contour = 10
    # Seleção de cada "subsample_factor"-ésimo ponto
    subsampled_lon = longitude[::subsample_factor]
    subsampled_lat = latitude[::subsample_factor]
    subsampled_u = vento_u_t0[::subsample_factor, ::subsample_factor]
    subsampled_v = vento_v_t0[::subsample_factor, ::subsample_factor]

    # Seleção de cada "subsample_factor_contour"-ésimo ponto
    subsampled_lon_contour = longitude[::subsample_factor_contour]
    subsampled_lat_contour = latitude[::subsample_factor_contour]
    subsampled_pressao_mmsl = pressao_mmsl_t0[::subsample_factor_contour, ::subsample_factor_contour]

    # Criar uma lista de níveis apenas com números pares
    pares = list(range(1012, 1000, -60))

     # Calculando a altura geopotencial para o tempo atual
    altura_geopotencial_t0 = geopot.isel(time=i)
    altura_em_metros = (altura_geopotencial_t0)/9.81

    # Restante do código permanece o mesmo
    contour = ax.contour(longitude, latitude, altura_em_metros, levels=10, colors='red', linestyles='dashed', linewidths=2)
  
    # Plotando as linhas de pressão ao nível médio do mar
    pressao_mmsl_t0 = pressao_mmsl_data.isel(time=i)
    contour_lines = ax.contour(subsampled_lon_contour, subsampled_lat_contour, subsampled_pressao_mmsl, levels=[1012], colors='k', linewidths=2, transform=ccrs.PlateCarree())
    contour_labels = ax.clabel(contour, inline=True, fontsize=10, fmt='%1.0f')
  
    # Adicionando etiquetas às linhas de contorno
    #plt.clabel(contour_lines, inline=True, fontsize=11, fmt='%1.0f hPa')

    img3 = ax.contour(longitude, latitude, pressao_mmsl_t0/100, colors='black', linewidths=1.2, levels=40)
    ax.clabel(img3, inline=1, inline_spacing=0, fontsize='10',fmt = '%1.0f', colors= 'black')
    plt.clabel(img3, inline=True, fontsize=10)

    # Adicionendo um círculo vermelho em Barra do Una
    ponto_lon = -45.7661
    ponto_lat = -23.7550
    ax.plot(ponto_lon, ponto_lat, 'ro', markersize=5, transform=ccrs.PlateCarree())

    # Criando uma paleta de cores personalizada
    azul_claro = 'lightskyblue'
    azul_escuro = 'darkblue'
    rosa_claro = 'lightpink'
    rosa_escuro = 'darkred'

    cmap_colors = [(0, azul_claro), (0.6, azul_escuro), (0.6, rosa_claro), (1, rosa_escuro)]
    custom_cmap = mcolors.LinearSegmentedColormap.from_list('custom', cmap_colors)

    # Plotando os dados com os limites atualizados e subamostrados
    contour_plot = ax.contourf(longitude, latitude, magnitude, transform=ccrs.PlateCarree(), cmap=custom_cmap, levels=np.arange(30, 110, 10))

    # Adicionando linhas de grade
    gl = ax.gridlines(draw_labels=True, dms=True)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlabels_bottom = False
    gl.ylabels_left = False

    # Adicionando uma barra de cores para o sombreamento
    cbar = plt.colorbar(contour_plot, ax=ax, orientation='vertical', pad=0.05, aspect=27, shrink=0.8)
    cbar.set_label('Velocidade do Vento (m s$^{{-1}}$)', fontsize=13)
    cbar.ax.set_position([0.76, 0.11, 0.1, 0.77])
    cbar.ax.tick_params(labelsize=12)

    # Plotando os dados
    data = vento250['time'].values[i]
    data_dt = datetime.utcfromtimestamp(data.tolist() / 1e9)  # Converter para objeto datetime
    data_str = data_dt.strftime('%Y%m%d%HZ')
    unidade_medida = "g kg$^{-1}$"
    titulo = f'Altura Geopotencial (m) em 500 hPa, Pressão ao Nível Médio do Mar (hPa) \ne Velocidade do Vento (m s$^{{-1}}$) em 250 hPa - {data_str}'
    ax.set_title(titulo, fontsize=16, position=(0.5, 3.5))
    ax.coastlines()
    ax.add_feature(cfeature.LAND, edgecolor='dimgray')
    ax.add_feature(cfeature.STATES)
    ax.set_xlim(longitude_min,longitude_max)
    ax.set_ylim(latitude_min,latitude_max)
    ax.set_xticks(np.arange(longitude_min, longitude_max, 10))
    ax.set_yticks(np.arange(latitude_min, latitude_max, 10))
    ax.set_xticklabels(np.arange(longitude_min, longitude_max, 10), fontsize=10)
    ax.set_yticklabels(np.arange(latitude_min, latitude_max, 10), fontsize=10)
    ax.yaxis.set_major_locator(MultipleLocator(10))

    # Salvar a imagem
    output = "local"
    fig.savefig(f'{output}/MSLP{data_str}_{i:03d}.png', bbox_inches='tight', pad_inches=0, dpi=300)
