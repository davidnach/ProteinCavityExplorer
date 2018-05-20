import { Component, OnInit, Input } from '@angular/core';
import * as d3 from 'd3';
export type Datum = {name: string, value: number};
@Component({
  selector: 'barcharts',
  templateUrl: './barcharts.component.html',
  styleUrls: ['./barcharts.component.css']
})
export class BarchartsComponent implements OnInit {
test: string;
static = 'here';
radius = 10;
height = 300;
width = 600;
/*data: Datum[] =  ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    .map((month: string) => ({
      name: month,
      value: Math.random() * 100
    }));*/
range = 100;
xScale: d3.ScaleBand<string> = null;
  yScale: d3.ScaleLinear<number, number> = null;
  transform = '';
  chartWidth = this.width;
  chartHeight = this.height;
  barHeights: number[] = [];
  barWidth = 0;
  xCoordinates: number[] = [];
	@Input() residueTypeCount : Map<string,number>;
  constructor() {
	this.test = '3';
	 }
	
  ngOnInit(){
	this.test = this.residueTypeCount.get('testing').toString();
  }
  ngOnChanges() { 
	this.test =  this.residueTypeCount.get('testing').toString();
 	console.log(this.residueTypeCount.get('testing'));
	console.log(this.test);
	this.chartHeight = this.height;
    this.chartWidth = this.width;
    /*this.xScale = d3.scaleBand()
      .domain(this.data.map((item: Datum)=>item.name)).range([0, this.chartWidth])
      .paddingInner(0.5);*/
    this.yScale = d3.scaleLinear()
      .domain([0, this.range])
      .range([this.chartHeight, 0]);

    this.barWidth = this.xScale.bandwidth();
    //this.barHeights = this.data.map((item: Datum) =>this.barHeight(item.value));
    this.barHeights = [1,2,3];
    //this.xCoordinates = this.data.map((item: Datum) => this.xScale(item.name));
    this.xCoordinates = [1,2,3];
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
