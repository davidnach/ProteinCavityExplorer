import { Component, OnInit, Input,ElementRef,ViewEncapsulation } from '@angular/core';
import { CommonModule } from "@angular/common";
import * as d3 from 'd3';
export type Datum = {name: string, value: number};
@Component({
  encapsulation: ViewEncapsulation.None,
  selector: 'barcharts',
  templateUrl: './barcharts.component.html',
  styleUrls: ['./barcharts.component.css']
})
export class BarchartsComponent  {
  data = [ ];
  svgExists: boolean = false;
  
  xScale: d3.ScaleBand<string> = null
  yScale: d3.ScaleLinear<number,number> = null
  transform = '';
  chartHeight = 300;
  aminoAcids = [];
  chartWidth = 600;
  //chartWidth = document.getElementById('bardisplay').offsetWidth;
  barHeights = [];
  barWidth = 0;
  xCoordinates = [];
  range = 100;
	@Input() residueTypeCount : Map<string,number>;
  constructor(private element: ElementRef) {
	
	 }
   ngOnChanges(){
    this.aminoAcids = [];
    var counts = [];
    for(let entry of Array.from(this.residueTypeCount.keys())) {
	this.aminoAcids.push(entry);
	counts.push(this.residueTypeCount.get(entry));
    }
    console.log(this.aminoAcids);
    this.xScale = d3.scaleBand()
      .domain(this.aminoAcids).range([0, this.chartWidth])
      .paddingInner(0.5);
    this.yScale = d3.scaleLinear()
      .domain([0, this.range])
      .range([this.chartHeight, 0]);

    this.barWidth = this.xScale.bandwidth();
    this.barHeights = counts.map(count => this.barHeight(count));
    console.log(this.barHeights); 
    this.xCoordinates = this.aminoAcids.map(aminoAcid => this.xScale(aminoAcid));
    console.log(this.xCoordinates);
    // use transform to flip the chart upside down, so the bars start from bottom
    this.transform = `scale(1, -1) translate(0, ${- this.chartHeight})`;
  }

  clampHeight(value: number) {
    if (value < 0) {
      return 0;
    }
    if (this.chartHeight <= 0) {
      return 0
    }
    if (value > this.chartHeight) {
      return this.chartHeight;
    }
    return value;
  }

  barHeight(value) {
    return this.clampHeight(this.chartHeight - this.yScale(value));
  }

}
