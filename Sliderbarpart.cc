#include "sliderbarpart.h"

void CreateImagePause(vtkSmartPointer<vtkImageData> image)
{
	// Specify the size of the image data
	image->SetDimensions(10, 10, 1);
	image->AllocateScalars(VTK_UNSIGNED_CHAR, 4);
	int* dims = image->GetDimensions();

	// Fill the image with
	for (int y = 0; y < dims[1]; y++)
	{
		for (int x = 0; x < dims[0]; x++)
		{
			unsigned char* pixel =
				static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			if (x < 4)
			{
				pixel[0] = 255;
				pixel[1] = 255;
				pixel[2] = 255;
				pixel[3] = 50;
			}
			else if (x > 5)
			{
				pixel[0] = 255;
				pixel[1] = 255;
				pixel[2] = 255;
				pixel[3] = 50;
			}
			else {
				pixel[3] = 0;
			}
		}
	}
}

void CreateImagePlay(vtkSmartPointer<vtkImageData> image)
{
	// Specify the size of the image data
	image->SetDimensions(200, 200, 1);
	image->AllocateScalars(VTK_UNSIGNED_CHAR, 4);
	int* dims = image->GetDimensions();

	// Fill the image with
	for (int x = 0; x < dims[0]; x++) {
		for (int y = 0; y< x / 2; y++) {
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			pixel[3] = 0;
		}
		for (int y = dims[1] - x / 2; y< dims[1]; y++) {
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			pixel[3] = 0;
		}
		for (int y = x / 2; y < dims[1] - x / 2; y++) {
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			pixel[0] = 255;
			pixel[1] = 255;
			pixel[2] = 255;
			pixel[3] = 50;
		}
	}
}

void CreateImageNextStep(vtkSmartPointer<vtkImageData> image)
{
	// Specify the size of the image data
	image->SetDimensions(200, 200, 1);
	image->AllocateScalars(VTK_UNSIGNED_CHAR, 4);
	int* dims = image->GetDimensions();

	for (int x = 0; x < 50; x++) {
		for (int y = 0; y < dims[1]; y++) {
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			pixel[0] = 255;
			pixel[1] = 255;
			pixel[2] = 255;
			pixel[3] = 50;
		}
	}
	for (int x = 50; x < 100; x++) {
		for (int y = 0; y < dims[1]; y++) {
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			pixel[3] = 0;
		}
	}
	// Fill the image with
	for (int x = 100; x < dims[0]; x++) {
		for (int y = 0; y < x - 100; y++) {
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			pixel[3] = 0;
		}
		for (int y = dims[1] + 100 - x; y< dims[1]; y++) {
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			pixel[3] = 0;
		}
		for (int y = x - 100; y < dims[1] + 100 - x; y++) {
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			pixel[0] = 255;
			pixel[1] = 255;
			pixel[2] = 255;
			pixel[3] = 50;
		}
	}
}

void CreateImagePrevStep(vtkSmartPointer<vtkImageData> image)
{
	// Specify the size of the image data
	image->SetDimensions(200, 200, 1);
	image->AllocateScalars(VTK_UNSIGNED_CHAR, 4);
	int* dims = image->GetDimensions();

	for (int x = 150; x < dims[0]; x++) {
		for (int y = 0; y < dims[1]; y++) {
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			pixel[0] = 255;
			pixel[1] = 255;
			pixel[2] = 255;
			pixel[3] = 50;
		}
	}

	for (int x = 100; x < 150; x++) {
		for (int y = 0; y < dims[1]; y++) {
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			pixel[3] = 0;
		}
	}
	// Fill the image with
	for (int x = 0; x < dims[0] / 2; x++) {
		for (int y = 0; y < -x + 100; y++) {
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			pixel[3] = 0;
		}
		for (int y = dims[1] - 100 + x; y< dims[1]; y++) {
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			pixel[3] = 0;
		}
		for (int y = -x + 100; y < dims[1] - 100 + x; y++) {
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			pixel[0] = 255;
			pixel[1] = 255;
			pixel[2] = 255;
			pixel[3] = 50;
		}
	}
}