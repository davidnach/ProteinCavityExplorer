import { Component, OnInit, Input,ElementRef,ViewEncapsulation } from '@angular/core';
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
  svg;
	@Input() residueTypeCount : Map<string,number>;
  constructor(private element: ElementRef) {
	
	 }
   ngOnChanges(){
	console.log(this.data);
	this.data = [];
	var current;
	console.log(this.residueTypeCount);
	//console.log(this.data);
	for(let entry of Array.from(this.residueTypeCount.keys())){
		current = {
			aminoAcid : entry,
			count : this.residueTypeCount.get(entry)
		}
		this.data.push(current);
	}
	console.log(this.data);	
     this.generateBarChart();
   }
   generateBarChart(){
	    console.log('here again');
            // set the dimensions and margins of the graph
            let margin = {top: 5, right: 20, bottom: 30, left: 40};
            let width = 600 - margin.left - margin.right;
            let height = 600 - margin.top - margin.bottom;

            //create svg
	    if(!this.svgExists) {
		
            	this.svg = d3.select(this.element.nativeElement).append('svg')
            	.attr('width', width + margin.left + margin.right)
            	.attr('height', height + margin.top + margin.bottom)
            	.style('background-color', '#efefef');
		this.svgExists;
	     } 
            //plot area

            let chart = this.svg.append("g")
            .attr('class', 'bar')
            .attr('transform', `translate(${margin.left}, ${margin.top})`);

            let xDomain = this.data.map(d => d.aminoAcid);
            let yDomain = [0, d3.max(this.data, d=> d.count)];

            // set the scale for data domain
            let x = d3.scaleBand()
                    .domain(xDomain)
                    .rangeRound([0, width])
		    .padding(0.2);
                    

            let y = d3.scaleLinear()
                    .domain(yDomain)
                    .range([height, 0]);

                    // add the x Axis
                    this.svg.append("g")
                        .attr('class', 'x axis')
                        .attr('transform', `translate(${margin.left}, ${margin.top + height})`)
                        .call(d3.axisBottom(x));

                    // add the y Axis
                    this.svg.append("g")
                        .attr('class', 'y axis')
                        .attr('transform', `translate(${margin.left}, ${margin.top})`)
                        .call(d3.axisLeft(y));

                    // plot chart with data
                    this.svg.selectAll("bar")
                        .data(this.data)
                        .enter().append("rect")
                        .attr("class", "bar")
                        .attr("x", function(d) { return margin.left + x(d.aminoAcid) ; })
                        .attr("width", x.bandwidth)
                        .attr("y", function(d) { return y(d.count); })
                        .attr("height", function(d) { return height - y(d.count); });
                }




}
