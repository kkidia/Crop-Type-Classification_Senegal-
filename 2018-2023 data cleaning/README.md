## 1. Data Digitization for 2023 


## 2. Data Cleaning for 2018, 2019, 2020 and 2023 years
### Completed Tasks:
- GEE Asset/Crop polygons
- Polygon area estimation for each crop type data.
- Planet NICFI Pixel count for each crops, each polygons
- Graphic representation for area vs Pixel count

### Next tasks
- Cloud masking for Planet (Algorithms)
- Image Normalization: Scale the pixel values (e.g., 0-255) to a range that is typical for neural network inputs, such as 0-1 or -1 to 1.
- Data Augmentation: To increase the diversity of the training data and prevent overfitting, apply transformations like rotations, translations, scaling, and horizontal flipping.
- Patch Extraction: For high-resolution images, it might be necessary to create smaller, manageable patches. This makes the training process more efficient and helps in handling large images during deployment.
- DL (FCNN) Trainings
- Parallel ML training

### [Access to analysis functions](https://drive.google.com/drive/folders/1-581wdLjY0_tf__913l5cY8tq4n4SqBV?usp=sharing)

1. [data_storage.py](https://drive.google.com/file/d/1-6_x0L6_yxaj3oxwmGJoYbn6luBgcnwX/view?usp=drive_link): - has all data assets and preprocessing 2023 data:
2. [process_raw_data.py](https://drive.google.com/file/d/1-9158gNZZzkJLlUvEiqUkq6S7cVNLMf4/view?usp=drive_link)- cleaning all 4 datasets to have same structure such as naming, Class, Subclass and Name.
3. [extract_s2_planet_ncfi_data.py ](to be uploaded)- has the class to process band values using Sentinel 2 and Planet NICFI basemaps

