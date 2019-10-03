# Multi Species Occupancy Models
Implement and exemplify the use of multi-species occupancy models (MSOM)


This project contains several builds of multi-species occupancy models and serves the purpose of

1. Demonstrating the how MSOM works 
2. Showing how standard ecological data can be organized to run the models  
3. Showing how flexible the model is for several types of ecological data 
4 - Creating a backbone for other users to improve the model (please create branches!)

# Structure

## Start

The file ExampleData1.csv contains a simple ecological dataframe that can be used as a model. The data is represented by a matrix with sites as rows and species as columns. The numbers filling the matrix represent species INCIDENCES (not abundance!). The number of incidences of a species in a site (number filling the matrix) could represent:

1 - number of pitfall traps in which a species was found in a site   
2 - number of observers that saw a bird in a site during a bird watching event   
3 - number of hours in a day in which a butterfly species was observed in a site   
4 - number of trees in a site where an epiphyte species was observed   
5 - number of nets in which a bat species was observed in a site   
6 - number of sections within a transect in which a plant or a termite species was observed   
7 - many other examples   

Note that in all cases, there is a cap on the maximum number of incidences in a site. For example, if a butterfly species was observed every hour during a 6 hour period, then the minimum observation of the species is zero and the maximum is 6. Therefore, the values filling the matrix MUST lay somewhere in the range from zero to 6. This does not mean that some species must have been observed in observed hours (i.e. the maximum value in the matrix can be lower than 6).

The last column of this data has the maximum number of incidences that could be obtained.

This data structure allows:
1 - Estimating the occupancy (true presence/absence) of a each species in each site   
2 - Estimating the detection (how easy it is to observe) of each species in each site   
3 - Estimating the variability among species in occupancy   
4 - Estimating the variability among species in detection   
5 - Estimating the number of undetected species   
6 - Estimating other community metrics (eg. similarity in species composition) when some species were not detected   
7 - Estimating how species occupancy and detection vary as a function of site covariate (if environmental covariates were measured in each site)   
8 - Estimating how species richness change with scale (from group of sites (region) to individual sites to individual replicates within sites).   

This data structure DOES NOT allows:

1 - Estimating how occupancy and detection of species respond to environmental co-variates measured WITHIN sites (eg. difference in forest cover from one pitfall to the next).


script Demonstration.R creates a simple ec



Markdown Cheatsheet
===================

- - - - 

# Heading 1 #

    Markup :  # Heading 1 #

    -OR-

    Markup :  ============= (below H1 text)

## Heading 2 ##

    Markup :  ## Heading 2 ##

    -OR-

    Markup: --------------- (below H2 text)

### Heading 3 ###

    Markup :  ### Heading 3 ###

#### Heading 4 ####

    Markup :  #### Heading 4 ####


Common text

    Markup :  Common text

_Emphasized text_

    Markup :  _Emphasized text_ or *Emphasized text*

~~Strikethrough text~~

    Markup :  ~~Strikethrough text~~

__Strong text__

    Markup :  __Strong text__ or **Strong text**

___Strong emphasized text___

    Markup :  ___Strong emphasized text___ or ***Strong emphasized text***

[Named Link](http://www.google.fr/ "Named link title") and http://www.google.fr/ or <http://example.com/>

    Markup :  [Named Link](http://www.google.fr/ "Named link title") and http://www.google.fr/ or <http://example.com/>

[heading-1](#heading-1 "Goto heading-1")
    
    Markup: [heading-1](#heading-1 "Goto heading-1")

Table, like this one :

First Header  | Second Header
------------- | -------------
Content Cell  | Content Cell
Content Cell  | Content Cell

```
First Header  | Second Header
------------- | -------------
Content Cell  | Content Cell
Content Cell  | Content Cell
```

`code()`

    Markup :  `code()`

```javascript
    var specificLanguage_code = 
    {
        "data": {
            "lookedUpPlatform": 1,
            "query": "Kasabian+Test+Transmission",
            "lookedUpItem": {
                "name": "Test Transmission",
                "artist": "Kasabian",
                "album": "Kasabian",
                "picture": null,
                "link": "http://open.spotify.com/track/5jhJur5n4fasblLSCOcrTp"
            }
        }
    }
```

    Markup : ```javascript
             ```

* Bullet list
    * Nested bullet
        * Sub-nested bullet etc
* Bullet list item 2

~~~
 Markup : * Bullet list
              * Nested bullet
                  * Sub-nested bullet etc
          * Bullet list item 2
~~~

1. A numbered list
    1. A nested numbered list
    2. Which is numbered
2. Which is numbered

~~~
 Markup : 1. A numbered list
              1. A nested numbered list
              2. Which is numbered
          2. Which is numbered
~~~

- [ ] An uncompleted task
- [x] A completed task

~~~
 Markup : - [ ] An uncompleted task
          - [x] A completed task
~~~

> Blockquote
>> Nested blockquote

    Markup :  > Blockquote
              >> Nested Blockquote

_Horizontal line :_
- - - -

    Markup :  - - - -

_Image with alt :_

![picture alt](http://www.brightlightpictures.com/assets/images/portfolio/thethaw_header.jpg "Title is optional")

    Markup : ![picture alt](http://www.brightlightpictures.com/assets/images/portfolio/thethaw_header.jpg "Title is optional")

Foldable text:

<details>
  <summary>Title 1</summary>
  <p>Content 1 Content 1 Content 1 Content 1 Content 1</p>
</details>
<details>
  <summary>Title 2</summary>
  <p>Content 2 Content 2 Content 2 Content 2 Content 2</p>
</details>

    Markup : <details>
               <summary>Title 1</summary>
               <p>Content 1 Content 1 Content 1 Content 1 Content 1</p>
             </details>

```html
<h3>HTML</h3>
<p> Some HTML code here </p>
```

Hotkey:

<kbd>⌘F</kbd>

<kbd>⇧⌘F</kbd>

    Markup : <kbd>⌘F</kbd>

Hotkey list:

| Key | Symbol |
| --- | --- |
| Option | ⌥ |
| Control | ⌃ |
| Command | ⌘ |
| Shift | ⇧ |
| Caps Lock | ⇪ |
| Tab | ⇥ |
| Esc | ⎋ |
| Power | ⌽ |
| Return | ↩ |
| Delete | ⌫ |
| Up | ↑ |
| Down | ↓ |
| Left | ← |
| Right | → |

Emoji:

:exclamation: Use emoji icons to enhance text. :+1:  Look up emoji codes at [emoji-cheat-sheet.com](http://emoji-cheat-sheet.com/)

    Markup : Code appears between colons :EMOJICODE:



