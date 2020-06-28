 # Chemistry definitions

The purpose of this folder is to give centralized, open source easy access to people who want to automatically pull information relevant to a sequencing protocol.

This folder contains all the data necessary for chemistry definitions.

## Adapters

The adapters folder contains fasta files used in the different protocols. They are also linked in the `definitions.json` file.

## Whitelists

The whitelists folder contains the list of barcodes used in the protocols defined in the `definitions.json` file.

## Definitions

The `definitions.json` file holds the information relevant to all the chemistries defined here.
The main purpose of the file is to allow fast, consistente access to similar information accross protocols from one source.

If you want to add your chemistry, please follow those guidelines:

Here is a template for a new entry:

```
    "shortname":{
        "Description": "",
        "adapters": [],
        "whitelist": "",
        "html": ""
    }
```

* shortname: is the key that is used to call the definition on the json file. The format is defined as follows: "contraction of compagny name" + "product name" + "version and/or specifics". Example: 10xchromiumFBV3 for "10x genomics" + "chromium" + "Feature barcoding V3"
* Description is a free field for a longer name of the protocol
* adapters is a list of paths pointing towards the adapters used in this protocol. The paths are relative to the folder containing the `definitions.json`
* whitelist is a path towards the whitelisted barcode list of the protocol. The path is relative to the folder containing the `definitions.json`
* html is a path towards the html document of the protocol. The path is relative to the folder containing the `definitions.json`

