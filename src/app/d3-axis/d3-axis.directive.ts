import {Directive,  Input,ElementRef } from '@angular/core';
import * as d3 from 'd3';
@Directive({
  selector: '[appD3Axis]'
})
export class D3AxisDirective {
  @Input() maxCount;
  @Input() scale: any;
  @Input() orientation: 'vertical' | 'horizontal' = 'horizontal';
  initialized = false;
  constructor(private el: ElementRef) {}

  drawAxis() {
    switch (this.orientation) {
      case 'horizontal':
        d3.select(this.el.nativeElement).call(d3.axisBottom(this.scale));
        break;
      case 'vertical':
        d3.select(this.el.nativeElement).call(d3.axisLeft(this.scale));
    }
  }

  ngAfterViewInit() {
    // all the Inputs will be set before this gets called.
    // D3 needs to wait for view init to modify it
    this.initialized = true;
    this.drawAxis();
  }

  ngOnChanges() {
    if (this.initialized) {
      this.drawAxis();
    }
  }
}
