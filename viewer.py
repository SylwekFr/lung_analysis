import matplotlib.pyplot as plt
import SimpleITK as sitk

def sitk_show(img, title=None, margin=0.05, dpi=40):
    """show slice in 2d https://pyscience.wordpress.com/2014/10/19/image-segmentation-with-python-and-simpleitk/"""
    nda = sitk.GetArrayFromImage(img)
    spacing = img.GetSpacing()
    figsize = (1 + margin) * nda.shape[0] / dpi, (1 + margin) * nda.shape[1] / dpi
    extent = (0, nda.shape[1] * spacing[1], nda.shape[0] * spacing[0], 0)
    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_axes([margin, margin, 1 - 2 * margin, 1 - 2 * margin])

    plt.set_cmap("gray")
    ax.imshow(nda, extent=extent, interpolation=None)

    if title:
        plt.title(title)

    # plt.show()
    #plt.savefig(str(title) + ".png")


def show_overlayed_images(img1, img2, title, folder, margin=0.05, dpi=40):
    plt.close('all')

    nda1 = sitk.GetArrayFromImage(img1)
    nda2 = sitk.GetArrayFromImage(img2)
    spacing = img1.GetSpacing()
    figsize = (1 + margin) * nda1.shape[1] / dpi, (1 + margin) * nda1.shape[0] / dpi
    # figsize = (1 + margin) * nda1.shape[0] / dpi, (1 + margin) * nda1.shape[1] / dpi
    extent = (0, nda1.shape[1] * spacing[1], nda1.shape[0] * spacing[0], 0)
    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_axes([margin, margin, 1 - 2 * margin, 1 - 2 * margin])
    plt.set_cmap("gray")
    ax.imshow(nda1, extent=extent, interpolation=None)
    ax.imshow(nda2, extent=extent, interpolation=None, alpha=0.5)
    if title:
        plt.title(title)
    # plt.show()
    #plt.savefig(folder + title + ".png")
    plt.close()






