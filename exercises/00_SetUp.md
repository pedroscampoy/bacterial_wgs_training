## Bacterial WGS training : Exercise 0

|**Title**| Working environment setup.|
|---------|-------------------------------------------|
|**Training dataset:**|  
|**Questions:**| <ul><li>How do I install the software for the course?</li><li>Where do I get the data for the exercises?</li></ul>|
|**Objectives**:|<ul><li>In this document we will cover the working environment setup for the exercises.</li></ul>|  
|**Time estimation**:| 5 min |
|**Key points**:|<ul><li>Each practical is designed to work on this folder structure, so make sure you follow these steps correctly.</li></ul>|

#### IMPORTANT: Make sure you understand and execute these commands in the right order.

Open a new terminal and navigate to your home directory if you are not already there:

```
pwd
cd
pwd
```

Create the project folder for the practises of the course:

```
mkdir -p Documents/wgs
```

Navigate to the directory:

```
cd Documents/wgs
```

Download git repository:

```
git clone https://github.com/BU-ISCIII/bacterial_wgs_training.git
```

Download training dataset:

```
wget https://github.com/BU-ISCIII/bacterial_wgs_training/releases/download/1.0/training_dataset_250k.tar.gz
tar -xvzf training_dataset_250k.tar.gz
rm -f training_dataset_250k.tar.gz
```

